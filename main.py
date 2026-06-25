import sys
import csv
import json
import gzip
import uuid
import yaml
import time
import zipfile
import threading
from pathlib import Path
from datetime import datetime
import streamlit as st
import click
import vcfpy
import myvariant

# --- Core NGS Pipeline Functions ---

from Bio import SeqIO
from Levenshtein import distance as edit_distance

def parse_fastq_basic(fastq_file):
    reads = bases = qsum = qcount = 0
    for record in SeqIO.parse(fastq_file, "fastq"):
        reads += 1
        bases += len(record.seq)
        quals = record.letter_annotations["phred_quality"]
        qsum += sum(quals)
        qcount += len(quals)
    avg_quality = qsum / qcount if qcount else 0
    return {
        'reads': reads,
        'total_bases': bases,
        'avg_quality': round(avg_quality, 2),
        'reads_per_mb': round(reads / (bases / 1e6), 0) if bases else 0,
    }

def parse_fasta(fasta_path):
    return {record.id: str(record.seq).upper() for record in SeqIO.parse(fasta_path, "fasta")}

def index_ref(ref_dict, k=11):
    index = {}
    for chrom, seq in ref_dict.items():
        for i in range(max(0, len(seq) - k + 1)):
            kmer = seq[i:i + k]
            if "N" not in kmer:
                index.setdefault(kmer, []).append((chrom, i))
    return index

REV_COMP = str.maketrans('ATGCATGCNn', 'TACGTACGNn')
def reverse_complement(seq):
    return seq.translate(REV_COMP)[::-1]

def candidate_positions(seq, ref_dict, ref_idx, k=11):
    if len(seq) < k:
        for chrom, ref_seq in ref_dict.items():
            for pos in range(len(ref_seq) - len(seq) + 1):
                yield chrom, pos
    else:
        yield from ref_idx.get(seq[:k], [])

def write_alignments_tsv(rows, out_path):
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["qname", "flag", "rname", "pos", "mapq", "cigar", "seq", "qual"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

def align_reads(fastq_file, reference_fasta, output_alignments=None, threads=4):
    ref_dict = parse_fasta(reference_fasta)
    ref_idx = index_ref(ref_dict, 11)
    alignments = []

    for record in SeqIO.parse(fastq_file, "fastq"):
        seq = str(record.seq).upper()
        qual = "".join(chr(q + 33) for q in record.letter_annotations["phred_quality"])
        best = None
        for flag, s in ((0, seq), (16, reverse_complement(seq))):
            for chrom, rpos in candidate_positions(s, ref_dict, ref_idx):
                ref_win = ref_dict[chrom][rpos:rpos + len(s)]
                if len(ref_win) == len(s):
                    dist = edit_distance(s, ref_win)
                    if best is None or dist < best[0]:
                        best = (dist, chrom, rpos + 1, flag)

        dist, chrom, pos, flag = best if best and best[0] <= 5 else (None, None, 0, 4)
        alignments.append({
            "qname": record.id,
            "flag": flag,
            "rname": chrom or "*",
            "pos": pos,
            "mapq": 60 if flag != 4 else 0,
            "cigar": f"{len(seq)}M" if flag != 4 else "*",
            "seq": seq,
            "qual": qual,
        })

    if output_alignments:
        write_alignments_tsv(alignments, output_alignments)
    return alignments

def call_variants(alignments, reference_fasta, output_vcf, min_depth=2):
    ref_dict = parse_fasta(reference_fasta)
    pileup = {}
    for r in alignments:
        if r["flag"] & 4:
            continue
        chrom = r["rname"]
        pos = r["pos"] - 1
        seq = r["seq"]

        for i in range(len(seq)):
            ref_pos = pos + i
            key = (chrom, ref_pos)
            if key not in pileup:
                pileup[key] = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
            base = seq[i].upper()
            if base in pileup[key]:
                pileup[key][base] += 1

    vcf_lines = [
        "##fileformat=VCFv4.2\n",
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    ]

    for (chrom, ref_pos), counts in sorted(pileup.items()):
        total = sum(counts.values())
        if total < min_depth:
            continue

        ref_base = ref_dict[chrom][ref_pos].upper()
        for base, count in counts.items():
            if base != ref_base and count >= min_depth and count / total >= 0.2:
                vcf_lines.append(f"{chrom}\t{ref_pos+1}\t.\t{ref_base}\t{base}\t60\tPASS\tDP={total}\n")

    with gzip.open(output_vcf, 'wt', encoding='utf-8') as f:
        f.writelines(vcf_lines)

MYVARIANT_FIELDS = ",".join([
    "_id",
    "snpeff.ann",
    "cadd.consequence",
    "cadd.gene.feature_id",
    "cadd.gene.gene_id",
    "cadd.gene.genename",
    "clinvar.gene.symbol",
    "clinvar.hgvs",
    "clinvar.rcv.clinical_significance",
    "dbnsfp.genename",
    "dbnsfp.ensembl.geneid",
    "dbnsfp.ensembl.transcriptid",
    "dbsnp.rsid",
    "cadd.phred",
    "gnomad_exome.af.af",
    "gnomad_genome.af.af",
])


def normalize_chromosome_for_myvariant(chrom):
    chrom = str(chrom)
    if chrom.lower().startswith("chr"):
        return chrom
    if chrom in {"MT", "M"}:
        return "chrM"
    return f"chr{chrom}"


def myvariant_hgvs_id(chrom, pos, ref, alt):
    return f"{normalize_chromosome_for_myvariant(chrom)}:g.{pos}{ref}>{alt}"


def first_scalar(value):
    if value in (None, ""):
        return None
    if isinstance(value, list):
        for item in value:
            scalar = first_scalar(item)
            if scalar not in (None, ""):
                return scalar
        return None
    if isinstance(value, dict):
        for item in value.values():
            scalar = first_scalar(item)
            if scalar not in (None, ""):
                return scalar
        return None
    return value


def nested_values(data, path):
    if data is None:
        return []
    if not path:
        return [data]
    key, rest = path[0], path[1:]
    if isinstance(data, list):
        values = []
        for item in data:
            values.extend(nested_values(item, path))
        return values
    if isinstance(data, dict) and key in data:
        return nested_values(data[key], rest)
    return []


def pick_annotation_value(data, *paths):
    for path in paths:
        for value in nested_values(data, path.split(".")):
            scalar = first_scalar(value)
            if scalar not in (None, ""):
                return scalar
    return None


def vcf_info_value(value):
    if value is None:
        return ""
    text = str(value)
    for old, new in ((";", ","), ("=", ":"), ("\t", "_"), ("\n", "_"), (" ", "_")):
        text = text.replace(old, new)
    return text


def build_myvariant_info(alt, annotation):
    if not annotation or annotation.get("notfound"):
        ann = f"{alt}|no_annotation_found|MODIFIER||||||||||||||"
        return {"ANN": ann, "MV_SOURCE": "myvariant.info", "MV_STATUS": "not_found"}

    gene = pick_annotation_value(annotation, "snpeff.ann.gene_name", "clinvar.gene.symbol", "cadd.gene.genename", "dbnsfp.genename") or ""
    gene_id = pick_annotation_value(annotation, "snpeff.ann.gene_id", "cadd.gene.gene_id", "dbnsfp.ensembl.geneid") or ""
    transcript = pick_annotation_value(annotation, "snpeff.ann.feature_id", "cadd.gene.feature_id", "dbnsfp.ensembl.transcriptid") or ""
    effect = pick_annotation_value(annotation, "snpeff.ann.effect", "cadd.consequence") or "sequence_variant"
    impact = pick_annotation_value(annotation, "snpeff.ann.putative_impact") or "MODIFIER"
    biotype = pick_annotation_value(annotation, "snpeff.ann.transcript_biotype") or ""
    rank = pick_annotation_value(annotation, "snpeff.ann.rank") or ""
    hgvs_c = pick_annotation_value(annotation, "snpeff.ann.hgvs_c", "clinvar.hgvs.coding") or ""
    hgvs_p = pick_annotation_value(annotation, "snpeff.ann.hgvs_p", "clinvar.hgvs.protein") or ""
    clinical_significance = pick_annotation_value(annotation, "clinvar.rcv.clinical_significance")
    cadd_phred = pick_annotation_value(annotation, "cadd.phred")
    allele_frequency = pick_annotation_value(annotation, "gnomad_genome.af.af", "gnomad_exome.af.af")

    ann = "|".join(vcf_info_value(value) for value in [
        alt,
        effect,
        impact,
        gene,
        gene_id,
        "transcript" if transcript else "sequence_feature",
        transcript,
        biotype,
        rank,
        hgvs_c,
        hgvs_p,
        "",
        "",
        "",
        "",
        "",
    ])

    info = {
        "ANN": ann,
        "MV_ID": annotation.get("_id", ""),
        "MV_SOURCE": "myvariant.info",
        "MV_STATUS": "found",
    }
    if gene:
        info["MV_GENE"] = gene
    if clinical_significance:
        info["MV_CLNSIG"] = clinical_significance
    if cadd_phred:
        info["MV_CADD_PHRED"] = cadd_phred
    if allele_frequency:
        info["MV_AF"] = allele_frequency
    return info


def annotate_variants(vcf_file, output_vcf, assembly="hg38"):
    variant_rows = []
    output_lines = []
    info_headers = [
        '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations from MyVariant.info, formatted as an ANN-compatible summary when SnpEff fields are available">\n',
        '##INFO=<ID=MV_ID,Number=1,Type=String,Description="MyVariant.info variant identifier">\n',
        '##INFO=<ID=MV_SOURCE,Number=1,Type=String,Description="Annotation source">\n',
        '##INFO=<ID=MV_STATUS,Number=1,Type=String,Description="MyVariant.info lookup status">\n',
        '##INFO=<ID=MV_GENE,Number=1,Type=String,Description="Gene symbol reported by MyVariant.info">\n',
        '##INFO=<ID=MV_CLNSIG,Number=1,Type=String,Description="ClinVar clinical significance reported by MyVariant.info">\n',
        '##INFO=<ID=MV_CADD_PHRED,Number=1,Type=String,Description="CADD PHRED score reported by MyVariant.info">\n',
        '##INFO=<ID=MV_AF,Number=1,Type=String,Description="gnomAD allele frequency reported by MyVariant.info">\n',
    ]
    inserted_headers = False

    with gzip.open(vcf_file, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#CHROM') and not inserted_headers:
                output_lines.extend(info_headers)
                inserted_headers = True
                output_lines.append(line)
                continue
            if line.startswith('#'):
                output_lines.append(line)
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 8:
                chrom, pos, _var_id, ref, alt, _qual, _filt, _info = fields[:8]
                query_id = myvariant_hgvs_id(chrom, pos, ref, alt)
                variant_rows.append((fields, query_id))
            else:
                output_lines.append(line)

    annotations = {}
    if variant_rows:
        mv = myvariant.MyVariantInfo()
        query_ids = [query_id for _fields, query_id in variant_rows]
        results = mv.getvariants(query_ids, fields=MYVARIANT_FIELDS, assembly=assembly)
        for query_id, result in zip(query_ids, results):
            annotations[query_id] = result

    for fields, query_id in variant_rows:
        alt = fields[4]
        info_values = build_myvariant_info(alt, annotations.get(query_id))
        new_info_parts = [fields[7]] if fields[7] not in ('.', '') else []
        new_info_parts.extend(
            f"{key}={vcf_info_value(value)}"
            for key, value in info_values.items()
            if value not in (None, "")
        )
        fields[7] = ";".join(new_info_parts) if new_info_parts else "."
        output_lines.append("\t".join(fields) + "\n")

    with gzip.open(output_vcf, 'wt', encoding='utf-8') as f:
        f.writelines(output_lines)

def parse_annotation(annotated_vcf):
    annotations = []
    for record in vcfpy.Reader.from_path(str(annotated_vcf)):
        annotation = record.INFO.get('ANN', '')
        if isinstance(annotation, list):
            annotation = annotation[0] if annotation else ''
        annotations.append({
            'chromosome': record.CHROM,
            'position': record.POS,
            'ref': record.REF,
            'alt': record.ALT[0].value if record.ALT else "",
            'quality': record.QUAL,
            'annotation': annotation,
            'gene': record.INFO.get('MV_GENE', ''),
            'clinical_significance': record.INFO.get('MV_CLNSIG', ''),
            'cadd_phred': record.INFO.get('MV_CADD_PHRED', ''),
            'allele_frequency': record.INFO.get('MV_AF', ''),
            'annotation_status': record.INFO.get('MV_STATUS', ''),
        })
    return annotations

def generate_markdown_report(output_dir, qc_stats, variant_count, annotation_data=None):
    summary = {
        'Total Reads': qc_stats.get('reads', 'N/A'),
        'Total Bases': qc_stats.get('total_bases', 'N/A'),
        'Avg Quality': qc_stats.get('avg_quality', 'N/A'),
        'Variants Called': variant_count,
        'Variants Annotated': len(annotation_data) if annotation_data else 0,
    }

    report = []
    report.append("# NGS Pipeline Analysis Report\n")
    report.append("## Summary Statistics\n")
    for k, v in summary.items():
        report.append(f"- **{k}**: {v}\n")

    report.append("\n## Quality Control Stats\n")
    report.append("```json\n" + json.dumps(qc_stats, indent=2) + "\n```\n")

    if annotation_data:
        report.append("\n## Variant Annotations (Top 10)\n")
        report.append("| Chromosome | Position | Ref | Alt | Gene | Effect | ClinVar | CADD PHRED | gnomAD AF | Status |\n")
        report.append("|------------|----------|-----|-----|------|--------|---------|------------|----------|--------|\n")
        for v in annotation_data[:10]:
            effect = v['annotation'].split('|')[1] if '|' in v['annotation'] else ''
            report.append(f"| {v['chromosome']} | {v['position']} | {v['ref']} | {v['alt']} | {v['gene']} | {effect} | {v['clinical_significance']} | {v['cadd_phred']} | {v['allele_frequency']} | {v['annotation_status']} |\n")

    report_path = Path(output_dir) / 'report.md'
    with open(report_path, 'w', encoding='utf-8') as f:
        f.writelines(report)
    return report_path

def run_full_pipeline(config, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fastq_input = config.get('fastq_input')
    reference = config.get('reference_genome')
    threads = config.get('threads', 4)
    enable_annotation = config.get('enable_annotation', False)
    annotation_assembly = config.get('annotation_assembly', 'hg38')

    # 1. QC
    qc_stats = parse_fastq_basic(fastq_input)

    # 2 & 3. Alignment
    alignments_file = output_dir / 'alignments.tsv'
    alignments = align_reads(fastq_input, reference, alignments_file, threads=threads)

    # 4. Variant Calling
    vcf_file = output_dir / 'variants.vcf.gz'
    call_variants(alignments, reference, vcf_file)

    # Count variants
    variant_count = 0
    with gzip.open(vcf_file, 'rt') as f:
        for line in f:
            if not line.startswith('#'):
                variant_count += 1

    # 5. Annotation
    annotation_data = None
    if enable_annotation:
        annotated_vcf = output_dir / 'variants_annotated.vcf.gz'
        annotate_variants(vcf_file, annotated_vcf, assembly=annotation_assembly)
        annotation_data = parse_annotation(annotated_vcf)

    # 6. Report
    generate_markdown_report(output_dir, qc_stats, variant_count, annotation_data)

# --- CLI Implementation ---

@click.group()
def cli():
    """NGS Data Processing Pipeline CLI"""
    pass

@cli.command()
@click.option('--config', required=True, type=click.Path(exists=True), help='Path to config YAML file')
@click.option('--output', required=True, type=click.Path(), help='Output directory')
def run(config, output):
    """Run the full NGS pipeline"""
    click.echo(f"Starting NGS Pipeline...")
    with open(config, 'r', encoding='utf-8') as f:
        config_data = yaml.safe_load(f)
    run_full_pipeline(config_data, output)
    click.echo(f"Pipeline completed successfully. Results saved to {output}")

# --- Streamlit Implementation ---

def run_streamlit_app():
    st.set_page_config(page_title="NGS Pipeline", page_icon="NGS", layout="wide")

    # Session state for navigation
    if "nav_selection" not in st.session_state:
        st.session_state["nav_selection"] = "Home"

    if "selected_job_view" not in st.session_state:
        st.session_state["selected_job_view"] = None

    RESULTS_DIR = Path("./web_results")
    RESULTS_DIR.mkdir(exist_ok=True)

    # Header Layout
    col_logo, col_nav = st.columns([2, 1])
    with col_logo:
        st.title("NGS Pipeline")
    with col_nav:
        nav_select = st.segmented_control("Navigation", options=["Home", "Submit", "Jobs"], default=st.session_state["nav_selection"])
        if nav_select != st.session_state["nav_selection"]:
            st.session_state["nav_selection"] = nav_select
            st.rerun()

    st.divider()

    # 1. HOME SCREEN
    if st.session_state["nav_selection"] == "Home":
        # Hero Section
        st.write("# NGS Pipeline")
        st.write("Enterprise bioinformatics platform for genomic sequencing analysis. Quality control, alignment, variant calling, and annotation in a single workflow.")

        col_act1, col_act2 = st.columns([1, 8])
        with col_act1:
            if st.button("Start Analysis", key="btn_start_analysis"):
                st.session_state["nav_selection"] = "Submit"
                st.rerun()
        with col_act2:
            st.button("Learn More", key="btn_learn_more", disabled=True)

        st.write("")
        st.write("")

        # Features Grid
        col_f1, col_f2, col_f3 = st.columns(3)
        with col_f1:
            with st.container(border=True):
                st.write("**PERFORMANCE**")
                st.write("### High Performance")
                st.write("Optimized algorithms and parallel processing for rapid analysis of large genomic datasets.")
        with col_f2:
            with st.container(border=True):
                st.write("**SCIENCE**")
                st.write("### Scientific Grade")
                st.write("Industry-standard tools and validated methodologies for research and clinical applications.")
        with col_f3:
            with st.container(border=True):
                st.write("**ANALYTICS**")
                st.write("### Detailed Reports")
                st.write("Comprehensive quality metrics, coverage analysis, and variant annotations in interactive reports.")

        st.write("")
        st.write("")

        # Stages Section
        st.write("## Analysis Pipeline")
        st.write("Six-stage workflow from raw sequencing data to annotated variants")

        col_s1, col_s2, col_s3, col_s4, col_s5, col_s6 = st.columns(6)
        stages_data = [
            ("01", "Quality Control", "Python Native"),
            ("02", "Index Setup", "In-Memory Index"),
            ("03", "Alignment", "Local Alignment"),
            ("04", "Variant Call", "Vectorized Pileup"),
            ("05", "Annotation", "Genomic Annotation"),
            ("06", "Report", "Markdown Output")
        ]
        columns_list = [col_s1, col_s2, col_s3, col_s4, col_s5, col_s6]
        for idx, (num, name, tool) in enumerate(stages_data):
            with columns_list[idx]:
                with st.container(border=True):
                    st.write(f"**{num}**")
                    st.write(f"**{name}**")
                    st.write(f"*{tool}*")

    # 2. SUBMIT SCREEN
    elif st.session_state["nav_selection"] == "Submit":
        st.write("## Submit New Analysis")
        st.write("Upload sequencing data and configure pipeline parameters")

        with st.container(border=True):
            fastq_file = st.file_uploader("Sequencing Reads (FASTQ)", type=["fastq", "fq"])
            ref_file = st.file_uploader("Reference Genome (FASTA)", type=["fasta", "fa"])
            threads = st.number_input("CPU Threads", min_value=1, max_value=128, value=4)
            enable_annotation = st.checkbox("Enable Variant Annotation (MyVariant.info)", value=True)
            annotation_assembly = st.selectbox("Annotation Assembly", options=["hg38", "hg19"], index=0)

            if st.button("Submit Analysis", key="btn_submit_job"):
                if not fastq_file or not ref_file:
                    st.error("Please select both a FASTQ file and a reference genome.")
                else:
                    job_id = str(uuid.uuid4())[:8]
                    job_dir = RESULTS_DIR / job_id
                    job_dir.mkdir(parents=True, exist_ok=True)

                    # Save files
                    fastq_path = job_dir / fastq_file.name
                    ref_path = job_dir / ref_file.name
                    with open(fastq_path, "wb") as f:
                        f.write(fastq_file.getbuffer())
                    with open(ref_path, "wb") as f:
                        f.write(ref_file.getbuffer())

                    # Config
                    config = {
                        'fastq_input': str(fastq_path),
                        'reference_genome': str(ref_path),
                        'threads': threads,
                        'enable_annotation': enable_annotation,
                        'annotation_assembly': annotation_assembly
                    }

                    # Start thread
                    t = threading.Thread(target=run_job_thread, args=(job_id, fastq_path, ref_path, config, RESULTS_DIR))
                    t.daemon = True
                    t.start()

                    st.success(f"Job submitted! ID: {job_id}")
                    time.sleep(1)
                    st.session_state["nav_selection"] = "Jobs"
                    st.rerun()

    # 3. JOBS SCREEN
    elif st.session_state["nav_selection"] == "Jobs":
        st.write("## Analysis History")
        st.write("View and manage your submitted analyses")

        job_dirs = list(RESULTS_DIR.glob('*'))
        jobs_data = []
        for d in job_dirs:
            info_file = d / "job_info.json"
            if info_file.exists():
                with open(info_file, "r", encoding='utf-8') as f:
                    meta = json.load(f)
                    jobs_data.append(meta)

        # Sort by creation date descending
        jobs_data = sorted(jobs_data, key=lambda x: x.get('created_at', ''), reverse=True)

        if not jobs_data:
            st.info("No analyses submitted yet.")
        else:
            for job in jobs_data:
                job_id = job['job_id']
                status = job['status']
                progress = job.get('progress', 0)
                created_at = job.get('created_at', '')

                with st.container(border=True):
                    col_info, col_status, col_prog, col_act = st.columns([2, 1, 2, 2])
                    with col_info:
                        st.write(f"**Analysis {job_id}**")
                        st.write(f"Submitted: {created_at}")
                        if status == 'failed' and 'error' in job:
                            st.write(f":red[Error: {job['error']}]")

                    with col_status:
                        if status == 'completed':
                            st.success("COMPLETED")
                        elif status == 'running':
                            st.info("RUNNING")
                        elif status == 'failed':
                            st.error("FAILED")
                        else:
                            st.warning("PENDING")

                    with col_prog:
                        st.progress(progress / 100)

                    with col_act:
                        if status == 'completed':
                            col_a1, col_a2 = st.columns(2)
                            with col_a1:
                                if st.button("View Report", key=f"view_{job_id}"):
                                    st.session_state["selected_job_view"] = job_id
                            with col_a2:
                                # Prepare ZIP
                                zip_path = RESULTS_DIR / f"{job_id}_results.zip"
                                if not zip_path.exists():
                                    job_dir = RESULTS_DIR / job_id
                                    with zipfile.ZipFile(zip_path, 'w') as zf:
                                        for f in job_dir.glob('*'):
                                            if f.is_file() and f.name != "job_info.json":
                                                zf.write(f, arcname=f.name)

                                if zip_path.exists():
                                    with open(zip_path, "rb") as f:
                                        st.download_button(
                                            label="Download ZIP",
                                            data=f.read(),
                                            file_name=f"{job_id}_results.zip",
                                            mime="application/zip",
                                            key=f"dl_{job_id}"
                                        )

            # Render selected report details
            if st.session_state["selected_job_view"]:
                view_id = st.session_state["selected_job_view"]
                st.write("---")
                st.write(f"### Report Details for Analysis {view_id}")
                report_path = RESULTS_DIR / view_id / "report.md"
                if report_path.exists():
                    with open(report_path, "r", encoding='utf-8') as f:
                        st.markdown(f.read())

    st.divider()
    st.write("NGS Pipeline v1")
    st.write("(c) 2026 Nikil Krishna")

def run_job_thread(job_id, fastq_path, ref_path, config, RESULTS_DIR):
    job_dir = RESULTS_DIR / job_id
    meta = {
        'job_id': job_id,
        'status': 'running',
        'progress': 10,
        'created_at': datetime.now().isoformat(),
        'config': {
            'threads': config['threads'],
            'enable_annotation': config['enable_annotation'],
            'annotation_assembly': config.get('annotation_assembly', 'hg38')
        }
    }
    try:
        with open(job_dir / "job_info.json", "w", encoding='utf-8') as f:
            json.dump(meta, f)

        # QC Stage
        meta['progress'] = 20
        with open(job_dir / "job_info.json", "w", encoding='utf-8') as f:
            json.dump(meta, f)
        qc_stats = parse_fastq_basic(fastq_path)

        # Index & Align
        meta['progress'] = 50
        with open(job_dir / "job_info.json", "w", encoding='utf-8') as f:
            json.dump(meta, f)
        alignments_file = job_dir / 'alignments.tsv'
        alignments = align_reads(fastq_path, ref_path, alignments_file, threads=config['threads'])

        # Variant Calling
        meta['progress'] = 70
        with open(job_dir / "job_info.json", "w", encoding='utf-8') as f:
            json.dump(meta, f)
        vcf_file = job_dir / 'variants.vcf.gz'
        call_variants(alignments, ref_path, vcf_file)

        # Count variants
        variant_count = 0
        with gzip.open(vcf_file, 'rt') as f:
            for line in f:
                if not line.startswith('#'):
                    variant_count += 1

        # Annotation
        meta['progress'] = 90
        with open(job_dir / "job_info.json", "w", encoding='utf-8') as f:
            json.dump(meta, f)
        annotation_data = None
        if config['enable_annotation']:
            annotated_vcf = job_dir / 'variants_annotated.vcf.gz'
            annotate_variants(vcf_file, annotated_vcf, assembly=config.get('annotation_assembly', 'hg38'))
            annotation_data = parse_annotation(annotated_vcf)

        # Report
        generate_markdown_report(job_dir, qc_stats, variant_count, annotation_data)
        meta['status'] = 'completed'
        meta['progress'] = 100
    except Exception as e:
        meta['status'] = 'failed'
        meta['error'] = str(e)
    finally:
        with open(job_dir / "job_info.json", "w", encoding='utf-8') as f:
            json.dump(meta, f)

# --- Launcher ---

if st.runtime.exists():
    run_streamlit_app()
else:
    if __name__ == '__main__':
        if len(sys.argv) > 1 and sys.argv[1] in ['run', '--help', '-h']:
            cli()
        else:
            app_path = str(Path(__file__).resolve())
            import streamlit.web.bootstrap as stbootstrap
            stbootstrap.run(app_path, False, [], flag_options={
                "server.headless": False,
                "browser.gatherUsageStats": False,
            })
