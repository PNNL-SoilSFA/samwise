#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * AuxModule_2_mvp.nf
 *
 * Runs selected MVP modules against assemblies from Module 2 and,
 * optionally, Module 2b.
 *
 * Example:
 *
 *   nextflow run AuxModule_2_mvp.nf \
 *       --working_dir /path/to/results \
 *       --mvp_modules "0,1,2,3,6" \
 *       --threads 24 \
 *       --memory_gb 192 \
 *       --genomad_db_path /path/to/databases/genomad_db \
 *       --checkv_db_path /path/to/databases/checkv_db \
 *       --install_databases false
 *
 * Supported MVP module selections:
 *
 *   0,1,2,3,4,5,6,7,99,100
 *
 * Modules are always executed in canonical numerical order, regardless of
 * the order supplied to --mvp_modules.
 *
 * MVP writes directly into:
 *
 *   <working_dir>/AuxModule_2_mvp/
 *
 * This preserves MVP's required directory names and allows completed
 * outputs to remain available if a later MVP module fails.
 */


/*
 * General inputs and output directories.
 */
params.working_dir = null
params.output_dir = null

params.assembly_manifest = null
params.trimmed_manifest = null

params.include_coassemblies = false
params.coassembly_manifest = null
params.coassembly_trimmed_manifest = null


/*
 * MVP module selection.
 *
 * Examples:
 *
 *   --mvp_modules "0,1,2,3,6"
 *   --mvp_modules "4,5,100"
 *   --mvp_modules "6,100"
 *
 * If earlier modules are omitted, their outputs must already exist in
 * params.outdir from a previous run.
 */
params.mvp_modules = "0,1,2,3,4,5,100"


/*
 * MVP installation.
 */
params.mvp_version = "1.1.4"
params.auto_install = true
params.tool_env_dir = null


/*
 * geNomad and CheckV databases.
 */
params.install_databases = false
params.mvp_database_dir = null
params.genomad_db_path = null
params.checkv_db_path = null


/*
 * Resources.
 */
params.threads = 4
params.memory_gb = 0


/*
 * MVP Module 00 options.
 */
params.skip_check_errors = false


/*
 * MVP Module 01 options.
 */
params.min_seq_size = 0
params.genomad_relaxed = false
params.genomad_conservative = false
params.skip_modify_headers = false


/*
 * MVP Module 02 options.
 */
params.viral_min_genes = 1
params.host_viral_genes_ratio = 1


/*
 * MVP Module 03 options.
 */
params.min_ani = 95
params.min_tcov = 85
params.min_qcov = 0
params.read_type = "short"
params.unfiltered_protein_file = false


/*
 * MVP Module 04 options.
 */
params.interleaved = false
params.delete_mapping_intermediates = true


/*
 * MVP Module 05 options.
 */
params.covered_fraction = null
params.normalization = "RPKM"
params.filtration = "conservative"


/*
 * MVP Module 06 options.
 */
params.functional_fasta_files = "representative"

params.phrogs_evalue = 0.01
params.phrogs_score = 60

params.pfam_evalue = 0.01
params.pfam_score = 50

params.functional_ads = true
params.ads_evalue = 0.01
params.ads_score = 60
params.ads_seqid = 30

params.functional_rdrp = false
params.rdrp_evalue = 0.01
params.rdrp_score = 50

params.functional_dram = true
params.delete_functional_intermediates = true


/*
 * MVP Module 07 options.
 */
params.binning_sample_group = null
params.read_mapping_sample_group = null
params.keep_bam = false
params.delete_binning_intermediates = false


/*
 * MVP Module 99 options.
 *
 * Required only when 99 is included in --mvp_modules.
 *
 * Valid steps:
 *   setup_metadata
 *   prep_submission
 *
 * A template is required for prep_submission.
 */
params.miuvig_identifier = null
params.miuvig_step = null
params.miuvig_template = null


/*
 * General behavior.
 */
params.force = false


/*
 * Derived directories.
 */
params.results_dir = params.working_dir
    ? params.working_dir
    : (params.output_dir ? params.output_dir : ".")

params.module1_outdir =
    "${params.results_dir}/module_1_readtrimming"

params.module2_outdir =
    "${params.results_dir}/module_2_readassembly"

params.module2b_outdir =
    "${params.results_dir}/module_2b_coassembly"

params.outdir =
    "${params.results_dir}/AuxModule_2_mvp"


/*
 * Convert a path to an absolute normalized path.
 */
def absPath(value) {
    return java.nio.file.Paths
        .get(value.toString())
        .toAbsolutePath()
        .normalize()
        .toString()
}


/*
 * Parse a comma-, semicolon-, or whitespace-separated module list.
 */
def parseMvpModules(value) {
    def aliases = [
        "00"  : "0",
        "01"  : "1",
        "02"  : "2",
        "03"  : "3",
        "04"  : "4",
        "05"  : "5",
        "06"  : "6",
        "07"  : "7",
        "099" : "99"
    ]

    return value
        .toString()
        .split(/[,;\s]+/)
        .collect { item -> item.trim() }
        .findAll { item -> item }
        .collect { item -> aliases.containsKey(item) ? aliases[item] : item }
        .unique()
}


workflow {

    def supported_modules = [
        "0",
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "99",
        "100"
    ]

    def canonical_order = [
        "0",
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "99",
        "100"
    ]

    def selected_modules =
        parseMvpModules(params.mvp_modules)

    if (!selected_modules) {
        error(
            """
            No MVP modules were selected.

            Use, for example:

              --mvp_modules "0,1,2,3,6"
            """.stripIndent()
        )
    }

    def invalid_modules =
        selected_modules.findAll { module ->
            !(module in supported_modules)
        }

    if (invalid_modules) {
        error(
            """
            Unsupported MVP module selection:

              ${invalid_modules.join(', ')}

            Supported modules:

              ${supported_modules.join(', ')}
            """.stripIndent()
        )
    }

    def ordered_modules =
        canonical_order.findAll { module ->
            module in selected_modules
        }

    if (!(params.read_type.toString() in ["short", "long"])) {
        error(
            """
            Invalid --read_type value: ${params.read_type}

            Supported values:
              short
              long
            """.stripIndent()
        )
    }

    if (!(params.normalization.toString() in ["RPKM", "FPKM"])) {
        error(
            """
            Invalid --normalization value: ${params.normalization}

            Supported values:
              RPKM
              FPKM
            """.stripIndent()
        )
    }

    if (
        !(
            params.filtration.toString() in
            ["relaxed", "conservative"]
        )
    ) {
        error(
            """
            Invalid --filtration value: ${params.filtration}

            Supported values:
              relaxed
              conservative
            """.stripIndent()
        )
    }

    if (
        params.genomad_relaxed.toString().toBoolean() &&
        params.genomad_conservative.toString().toBoolean()
    ) {
        error(
            """
            --genomad_relaxed and --genomad_conservative cannot
            both be enabled.
            """.stripIndent()
        )
    }

    if ("99" in selected_modules) {

        if (!params.miuvig_identifier) {
            error(
                """
                MVP Module 99 was selected, but no MIUViG identifier
                was supplied.

                Provide:

                  --miuvig_identifier IDENTIFIER
                """.stripIndent()
            )
        }

        if (
            !params.miuvig_step ||
            !(
                params.miuvig_step.toString() in
                ["setup_metadata", "prep_submission"]
            )
        ) {
            error(
                """
                MVP Module 99 requires:

                  --miuvig_step setup_metadata

                or:

                  --miuvig_step prep_submission
                """.stripIndent()
            )
        }

        if (
            params.miuvig_step.toString() == "prep_submission" &&
            !params.miuvig_template
        ) {
            error(
                """
                MVP Module 99 prep_submission requires:

                  --miuvig_template /path/to/template.tsv
                """.stripIndent()
            )
        }
    }

    def assembly_manifest_file =
        params.assembly_manifest ?:
        "${params.module2_outdir}/summary/assembly_manifest.tsv"

    def trimmed_manifest_file =
        params.trimmed_manifest ?:
        "${params.module1_outdir}/summary/trimmed_manifest.tsv"

    def include_coassemblies =
        params.include_coassemblies
            .toString()
            .toBoolean()

    log.info(
        "Auxiliary Module 2 MVP output directory: " +
        "${params.outdir}"
    )

    log.info(
        "Selected MVP modules: " +
        "${ordered_modules.join(',')}"
    )

    log.info(
        "MVP modules will execute in this order: " +
        "${ordered_modules.join(' -> ')}"
    )

    log.info(
        "Assembly manifest: ${assembly_manifest_file}"
    )

    log.info(
        "Trimmed-read manifest: ${trimmed_manifest_file}"
    )

    log.info(
        "Include coassemblies: ${include_coassemblies}"
    )

    log.info(
        "MVP version: ${params.mvp_version}"
    )

    log.info(
        "MVP threads: ${params.threads}"
    )

    log.info(
        "MVP read type: ${params.read_type}"
    )

    if (!("0" in selected_modules)) {
        log.warn(
            "MVP Module 00 was not selected. Its required setup " +
            "and directory outputs must already exist."
        )
    }

    if (
        ("2" in selected_modules) &&
        !("1" in selected_modules)
    ) {
        log.warn(
            "MVP Module 02 was selected without Module 01. " +
            "Module 01 outputs must already exist."
        )
    }

    if (
        ("3" in selected_modules) &&
        !("2" in selected_modules)
    ) {
        log.warn(
            "MVP Module 03 was selected without Module 02. " +
            "Module 02 outputs must already exist."
        )
    }

    if (
        ("4" in selected_modules) &&
        !("3" in selected_modules)
    ) {
        log.warn(
            "MVP Module 04 was selected without Module 03. " +
            "Module 03 outputs and the mapping index must already exist."
        )
    }

    if (
        ("5" in selected_modules) &&
        !("4" in selected_modules)
    ) {
        log.warn(
            "MVP Module 05 was selected without Module 04. " +
            "Module 04 coverage outputs must already exist."
        )
    }

    if (
        ("6" in selected_modules) &&
        !("3" in selected_modules)
    ) {
        log.warn(
            "MVP Module 06 was selected without Module 03. " +
            "Module 03 protein and sequence outputs must already exist."
        )
    }

    if (
        ("7" in selected_modules) &&
        !("3" in selected_modules)
    ) {
        log.warn(
            "MVP Module 07 was selected without Module 03. " +
            "Module 03 outputs must already exist."
        )
    }

    def assembly_manifest_ch = channel.fromPath(
        assembly_manifest_file,
        type: "file",
        checkIfExists: true
    )

    def trimmed_manifest_ch = channel.fromPath(
        trimmed_manifest_file,
        type: "file",
        checkIfExists: true
    )

    def coassembly_manifest_ch
    def coassembly_reads_ch

    if (include_coassemblies) {

        def coassembly_manifest_file =
            params.coassembly_manifest ?:
            "${params.module2b_outdir}/summary/assembly_manifest.tsv"

        def coassembly_reads_file =
            params.coassembly_trimmed_manifest ?:
            "${params.module2b_outdir}/summary/coassembly_trimmed_manifest.tsv"

        log.info(
            "Coassembly assembly manifest: " +
            "${coassembly_manifest_file}"
        )

        log.info(
            "Coassembly trimmed-read manifest: " +
            "${coassembly_reads_file}"
        )

        coassembly_manifest_ch = channel.fromPath(
            coassembly_manifest_file,
            type: "file",
            checkIfExists: true
        )

        coassembly_reads_ch = channel.fromPath(
            coassembly_reads_file,
            type: "file",
            checkIfExists: true
        )
    }
    else {

        CREATE_EMPTY_COASSEMBLY_MANIFESTS()

        coassembly_manifest_ch =
            CREATE_EMPTY_COASSEMBLY_MANIFESTS
                .out
                .assembly_manifest

        coassembly_reads_ch =
            CREATE_EMPTY_COASSEMBLY_MANIFESTS
                .out
                .trimmed_manifest
    }

    SETUP_AUXMODULE2_MVP()

    PREPARE_MVP_METADATA(
        assembly_manifest_ch,
        trimmed_manifest_ch,
        coassembly_manifest_ch,
        coassembly_reads_ch
    )

    RUN_MVP(
        PREPARE_MVP_METADATA.out.metadata,
        PREPARE_MVP_METADATA.out.input_summary,
        SETUP_AUXMODULE2_MVP.out.status
    )
}


process CREATE_EMPTY_COASSEMBLY_MANIFESTS {

    tag "create_empty_coassembly_manifests"

    output:

    path "empty_coassembly_assembly_manifest.tsv",
        emit: assembly_manifest

    path "empty_coassembly_trimmed_manifest.tsv",
        emit: trimmed_manifest

    script:

    """
    set -euo pipefail

    printf '%s\\n' \
      'sample_id\tsafe_sample_id\tassembly_sample_id\tassembler\tassembly_mode\trarefaction_label\tassembly_strategy\trenamed_fasta' \
      > empty_coassembly_assembly_manifest.tsv

    printf '%s\\n' \
      'sample_id\tsafe_sample_id\tlayout\tread1\tread2\tinterleaved\tmerged\tfastp_html\tfastp_json' \
      > empty_coassembly_trimmed_manifest.tsv
    """
}


process SETUP_AUXMODULE2_MVP {

    tag "setup_auxmodule2_mvp"

    publishDir "${params.outdir}/setup",
        mode: "copy",
        pattern: "auxmodule2_mvp_tools_status.env"

    output:

    path "auxmodule2_mvp_tools_status.env",
        emit: status

    script:

    def env_dir = params.tool_env_dir
        ? absPath(params.tool_env_dir)
        : absPath("${params.outdir}/conda_envs/mvp")

    """
    set -euo pipefail

    STATUS_FILE="auxmodule2_mvp_tools_status.env"
    TOOL_ENV="${env_dir}"

    echo "Auxiliary Module 2 MVP setup started: \$(date)" \
        > "\$STATUS_FILE"

    echo "Requested MVP version: ${params.mvp_version}" \
        >> "\$STATUS_FILE"

    echo "Requested tool environment: \$TOOL_ENV" \
        >> "\$STATUS_FILE"

    echo "----------------------------------------" \
        >> "\$STATUS_FILE"

    check_mvp() {
        local prefix="\$1"

        if [[ "\$prefix" == "SYSTEM" ]]; then

            if ! command -v mvip >/dev/null 2>&1; then
                return 1
            fi

            mvip -h >> "\$STATUS_FILE" 2>&1 || return 1
            return 0
        fi

        if [[ ! -x "\$prefix/bin/mvip" ]]; then
            return 1
        fi

        "\$prefix/bin/mvip" -h \
            >> "\$STATUS_FILE" 2>&1 || return 1

        return 0
    }

    if [[ -d "\$TOOL_ENV" ]]; then

        echo "Existing MVP environment detected: \$TOOL_ENV" \
            >> "\$STATUS_FILE"

        if check_mvp "\$TOOL_ENV"; then

            echo "Existing MVP environment passed validation." \
                >> "\$STATUS_FILE"

            echo "TOOL_ENV=\$TOOL_ENV" \
                >> "\$STATUS_FILE"

            echo "Auxiliary Module 2 MVP setup finished: \$(date)" \
                >> "\$STATUS_FILE"

            exit 0
        fi

        echo "Existing MVP environment is incomplete or broken." \
            >> "\$STATUS_FILE"

        echo "Removing incomplete MVP environment." \
            >> "\$STATUS_FILE"

        rm -rf "\$TOOL_ENV"
    fi

    echo "Checking system/runtime MVP installation..." \
        >> "\$STATUS_FILE"

    if check_mvp "SYSTEM"; then

        echo "MVP is available from the system/runtime PATH." \
            >> "\$STATUS_FILE"

        echo "TOOL_ENV=SYSTEM" \
            >> "\$STATUS_FILE"

        echo "Auxiliary Module 2 MVP setup finished: \$(date)" \
            >> "\$STATUS_FILE"

        exit 0
    fi

    echo "MVP is not available from the system/runtime PATH." \
        >> "\$STATUS_FILE"

    if [[ "${params.auto_install}" != "true" ]]; then

        echo "ERROR: MVP is unavailable and auto_install is false." \
            >> "\$STATUS_FILE"

        exit 1
    fi

    INSTALLER=""

    if command -v mamba >/dev/null 2>&1; then

        INSTALLER="mamba"

        echo "Using mamba: \$(command -v mamba)" \
            >> "\$STATUS_FILE"

    elif command -v conda >/dev/null 2>&1; then

        INSTALLER="conda"

        echo "Using conda: \$(command -v conda)" \
            >> "\$STATUS_FILE"

    else

        echo "ERROR: Neither mamba nor conda was found in PATH." \
            >> "\$STATUS_FILE"

        exit 1
    fi

    mkdir -p "\$(dirname "\$TOOL_ENV")"

    echo "Creating MVP environment: \$TOOL_ENV" \
        >> "\$STATUS_FILE"

    "\$INSTALLER" create -y \
        -p "\$TOOL_ENV" \
        -c conda-forge \
        -c bioconda \
        "mvip=${params.mvp_version}" \
        >> "\$STATUS_FILE" 2>&1

    if ! check_mvp "\$TOOL_ENV"; then

        echo "ERROR: Newly created MVP environment failed validation." \
            >> "\$STATUS_FILE"

        exit 1
    fi

    echo "TOOL_ENV=\$TOOL_ENV" \
        >> "\$STATUS_FILE"

    echo "Auxiliary Module 2 MVP setup finished: \$(date)" \
        >> "\$STATUS_FILE"
    """
}


process PREPARE_MVP_METADATA {

    tag "prepare_mvp_metadata"

    publishDir "${params.outdir}/metadata",
        mode: "copy",
        pattern: "mvp_metadata.tsv"

    publishDir "${params.outdir}/metadata",
        mode: "copy",
        pattern: "mvp_input_summary.tsv"

    publishDir "${params.outdir}/logs",
        mode: "copy",
        pattern: "prepare_mvp_metadata.log"

    input:

    path assembly_manifest
    path trimmed_manifest
    path coassembly_manifest
    path coassembly_trimmed_manifest

    output:

    path "mvp_metadata.tsv",
        emit: metadata

    path "mvp_input_summary.tsv",
        emit: input_summary

    path "prepare_mvp_metadata.log",
        emit: log_file

    script:

    """
    set -euo pipefail

    python3 - \
        "${assembly_manifest}" \
        "${trimmed_manifest}" \
        "${coassembly_manifest}" \
        "${coassembly_trimmed_manifest}" \
        "${workflow.launchDir}" \
        "mvp_metadata.tsv" \
        "mvp_input_summary.tsv" \
        "prepare_mvp_metadata.log" <<'PY'
import csv
import re
import sys
from pathlib import Path


(
    assembly_manifest,
    trimmed_manifest,
    coassembly_manifest,
    coassembly_trimmed_manifest,
    launch_dir,
    output_metadata,
    output_summary,
    log_file,
) = sys.argv[1:]


launch_dir = Path(launch_dir).resolve()
output_metadata = Path(output_metadata)
output_summary = Path(output_summary)
log_file = Path(log_file)


def log(message):
    with log_file.open("a") as handle:
        print(message, file=handle)


def resolve_path(value):
    value = str(value or "").strip()

    if not value:
        return ""

    path = Path(value)

    if not path.is_absolute():
        path = launch_dir / path

    return str(path.resolve())


def safe_name(value):
    value = str(value or "").strip()
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    value = value.strip("._-")

    return value or "unnamed"


def read_tsv(path, required_columns):
    path = Path(path)

    if not path.exists():
        raise RuntimeError(
            f"Input manifest does not exist: {path}"
        )

    with path.open() as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\\t",
        )

        fields = set(reader.fieldnames or [])
        missing = set(required_columns) - fields

        if missing:
            raise RuntimeError(
                f"Manifest {path} is missing required columns: "
                f"{sorted(missing)}"
            )

        return list(reader)


def build_read_lookup(rows, label):
    lookup = {}

    for row in rows:
        sample_id = row["sample_id"].strip()
        safe_sample_id = row["safe_sample_id"].strip()
        layout = row["layout"].strip()

        if not sample_id:
            raise RuntimeError(
                f"{label} contains an empty sample_id"
            )

        if layout == "paired":
            read_path = resolve_path(row.get("read1", ""))
            read2_path = resolve_path(row.get("read2", ""))

            if not read_path or not read2_path:
                raise RuntimeError(
                    f"Paired sample {sample_id} is missing "
                    f"read1 or read2"
                )

            if not Path(read_path).is_file():
                raise RuntimeError(
                    f"Read 1 does not exist for sample "
                    f"{sample_id}: {read_path}"
                )

            if not Path(read2_path).is_file():
                raise RuntimeError(
                    f"Read 2 does not exist for sample "
                    f"{sample_id}: {read2_path}"
                )

        elif layout == "interleaved":
            read_path = resolve_path(
                row.get("interleaved", "")
            )

            if not read_path:
                raise RuntimeError(
                    f"Interleaved sample {sample_id} has no "
                    f"read path"
                )

            if not Path(read_path).is_file():
                raise RuntimeError(
                    f"Interleaved reads do not exist for sample "
                    f"{sample_id}: {read_path}"
                )

        else:
            raise RuntimeError(
                f"Unsupported read layout for sample "
                f"{sample_id}: {layout}"
            )

        record = {
            "sample_id": sample_id,
            "safe_sample_id": safe_sample_id,
            "layout": layout,
            "read_path": read_path,
        }

        for key in {sample_id, safe_sample_id}:

            if not key:
                continue

            if key in lookup:
                previous = lookup[key]

                if previous["sample_id"] != sample_id:
                    raise RuntimeError(
                        f"Ambiguous read lookup key in "
                        f"{label}: {key}"
                    )

            lookup[key] = record

    return lookup


log("Preparing MVP metadata")
log(f"Launch directory: {launch_dir}")
log(f"Assembly manifest: {assembly_manifest}")
log(f"Trimmed-read manifest: {trimmed_manifest}")
log(f"Coassembly manifest: {coassembly_manifest}")
log(
    "Coassembly trimmed-read manifest: "
    f"{coassembly_trimmed_manifest}"
)


required_read_columns = [
    "sample_id",
    "safe_sample_id",
    "layout",
    "read1",
    "read2",
    "interleaved",
]


trimmed_rows = read_tsv(
    trimmed_manifest,
    required_read_columns,
)

coassembly_read_rows = read_tsv(
    coassembly_trimmed_manifest,
    required_read_columns,
)


read_lookup = build_read_lookup(
    trimmed_rows,
    "trimmed-read manifest",
)

coassembly_read_lookup = build_read_lookup(
    coassembly_read_rows,
    "coassembly trimmed-read manifest",
)


required_assembly_columns = [
    "sample_id",
    "safe_sample_id",
    "assembly_sample_id",
    "assembler",
    "assembly_mode",
    "rarefaction_label",
    "assembly_strategy",
    "renamed_fasta",
]


ordinary_assemblies = read_tsv(
    assembly_manifest,
    required_assembly_columns,
)

coassemblies = read_tsv(
    coassembly_manifest,
    required_assembly_columns,
)


records = []
used_mvp_names = set()
skipped_assemblies = []


def add_assembly(row, source, lookup):
    sample_id = row["sample_id"].strip()
    safe_sample_id = row["safe_sample_id"].strip()
    assembly_sample_id = row["assembly_sample_id"].strip()
    assembler = row["assembler"].strip()
    assembly_mode = row["assembly_mode"].strip()
    rarefaction_label = row["rarefaction_label"].strip()
    assembly_strategy = row["assembly_strategy"].strip()
    assembly_path = resolve_path(row["renamed_fasta"])

    if not sample_id:
        raise RuntimeError(
            f"An assembly record from {source} has no sample_id"
        )

    if not assembly_path:
        raise RuntimeError(
            f"Missing assembly path for {source} sample "
            f"{sample_id}"
        )

    assembly_file = Path(assembly_path)

    if not assembly_file.is_file():
        raise RuntimeError(
            f"Assembly file does not exist for {sample_id}: "
            f"{assembly_path}"
        )

    if assembly_file.stat().st_size == 0:
        skipped_assemblies.append(
            (
                sample_id,
                assembler,
                assembly_mode,
                assembly_path,
                "empty assembly file",
            )
        )

        log(
            f"WARNING: skipping empty assembly for "
            f"{sample_id}: {assembly_path}"
        )

        return

    read_record = (
        lookup.get(sample_id) or
        lookup.get(safe_sample_id)
    )

    if read_record is None:
        raise RuntimeError(
            f"No trimmed-read record matched {source} "
            f"assembly {sample_id!r}/{safe_sample_id!r}"
        )

    components = [
        safe_sample_id or sample_id,
        assembler or "assembler",
        assembly_mode or "assembly",
    ]

    if rarefaction_label:
        components.append(rarefaction_label)

    if assembly_strategy:
        components.append(assembly_strategy)

    base_name = safe_name("_".join(components))
    mvp_name = base_name
    suffix = 1

    while mvp_name in used_mvp_names:
        suffix += 1
        mvp_name = f"{base_name}_{suffix}"

    used_mvp_names.add(mvp_name)

    records.append(
        {
            "Sample_number": len(records) + 1,
            "Sample": mvp_name,
            "Assembly_Path": assembly_path,
            "Read_Path": read_record["read_path"],
            "Variable": source,
            "source_sample_id": sample_id,
            "safe_sample_id": safe_sample_id,
            "assembly_sample_id": assembly_sample_id,
            "assembler": assembler,
            "assembly_mode": assembly_mode,
            "rarefaction_label": rarefaction_label,
            "assembly_strategy": assembly_strategy,
            "read_layout": read_record["layout"],
        }
    )


for row in ordinary_assemblies:
    add_assembly(
        row=row,
        source="individual",
        lookup=read_lookup,
    )


for row in coassemblies:
    add_assembly(
        row=row,
        source="coassembly",
        lookup=coassembly_read_lookup,
    )


if not records:
    raise RuntimeError(
        "No usable MVP input records were generated"
    )


with output_metadata.open("w") as out:
    fields = [
        "Sample_number",
        "Sample",
        "Assembly_Path",
        "Read_Path",
        "Variable",
    ]

    writer = csv.DictWriter(
        out,
        delimiter="\\t",
        lineterminator="\\n",
        fieldnames=fields,
    )

    writer.writeheader()

    for record in records:
        writer.writerow(
            {
                key: record[key]
                for key in fields
            }
        )


with output_summary.open("w") as out:
    fields = [
        "Sample_number",
        "Sample",
        "source_sample_id",
        "safe_sample_id",
        "assembly_sample_id",
        "assembler",
        "assembly_mode",
        "rarefaction_label",
        "assembly_strategy",
        "read_layout",
        "Assembly_Path",
        "Read_Path",
        "Variable",
    ]

    writer = csv.DictWriter(
        out,
        delimiter="\\t",
        lineterminator="\\n",
        fieldnames=fields,
    )

    writer.writeheader()
    writer.writerows(records)


log(f"MVP input records written: {len(records)}")

log(
    "Individual assembly records: "
    f"{sum(r['Variable'] == 'individual' for r in records)}"
)

log(
    "Coassembly records: "
    f"{sum(r['Variable'] == 'coassembly' for r in records)}"
)

log(
    "Empty assembly records skipped: "
    f"{len(skipped_assemblies)}"
)

for skipped in skipped_assemblies:
    log(
        "Skipped assembly: "
        + "\\t".join(str(value) for value in skipped)
    )
PY
    """
}


process RUN_MVP {

    tag "run_auxmodule2_mvp"

    stageInMode "symlink"

    /*
     * MVP writes directly into params.outdir.
     *
     * Do not add mvp_work as a publishDir output. The direct-write
     * behavior preserves MVP's required directory structure and retains
     * partial results if a later module fails.
     */
    publishDir "${params.outdir}/logs",
        mode: "copy",
        pattern: "mvp_commands.log"

    publishDir "${params.outdir}/logs",
        mode: "copy",
        pattern: "mvp_run.log"

    publishDir "${params.outdir}/summary",
        mode: "copy",
        pattern: "mvp_input_summary.tsv"

    publishDir "${params.outdir}/summary",
        mode: "copy",
        pattern: "mvp_complete.txt"

    cpus {
        params.threads as int
    }

    memory {
        def gb = params.memory_gb as int

        return gb > 0
            ? "${gb} GB"
            : null
    }

    input:

    path metadata
    path input_summary
    path tools_status

    output:

    path "mvp_commands.log",
        emit: commands

    path "mvp_run.log",
        emit: log_file

    path "mvp_input_summary.tsv",
        emit: input_summary

    path "mvp_complete.txt",
        emit: completion_marker

    script:

    def selected_modules =
        parseMvpModules(params.mvp_modules)

    def canonical_order = [
        "0",
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "99",
        "100"
    ]

    def ordered_modules =
        canonical_order.findAll { module ->
            module in selected_modules
        }

    def run_module_0 =
        "0" in selected_modules

    def run_module_1 =
        "1" in selected_modules

    def run_module_2 =
        "2" in selected_modules

    def run_module_3 =
        "3" in selected_modules

    def run_module_4 =
        "4" in selected_modules

    def run_module_5 =
        "5" in selected_modules

    def run_module_6 =
        "6" in selected_modules

    def run_module_7 =
        "7" in selected_modules

    def run_module_99 =
        "99" in selected_modules

    def run_module_100 =
        "100" in selected_modules

    def install_databases =
        params.install_databases
            .toString()
            .toBoolean()

    def force =
        params.force
            .toString()
            .toBoolean()

    def mvp_work =
        absPath(params.outdir)

    def database_root = params.mvp_database_dir
        ? absPath(params.mvp_database_dir)
        : "${mvp_work}/databases"

    def genomad_database = params.genomad_db_path
        ? absPath(params.genomad_db_path)
        : "${database_root}/genomad_db"

    def checkv_database = params.checkv_db_path
        ? absPath(params.checkv_db_path)
        : "${database_root}/checkv_db"

    def module00_database_flag =
        install_databases
            ? ""
            : "--skip_install_databases"

    def skip_check_errors_argument =
        params.skip_check_errors
            .toString()
            .toBoolean()
            ? "--skip_check_errors"
            : ""

    def min_seq_argument =
        (params.min_seq_size as int) > 0
            ? "--min_seq_size ${params.min_seq_size}"
            : ""

    def genomad_filter_argument = ""

    if (
        params.genomad_relaxed
            .toString()
            .toBoolean()
    ) {
        genomad_filter_argument =
            "--genomad_relaxed"
    }
    else if (
        params.genomad_conservative
            .toString()
            .toBoolean()
    ) {
        genomad_filter_argument =
            "--genomad_conservative"
    }

    def skip_modify_headers_argument =
        params.skip_modify_headers
            .toString()
            .toBoolean()
            ? "--skip_modify_headers"
            : ""

    def unfiltered_protein_argument =
        params.unfiltered_protein_file
            .toString()
            .toBoolean()
            ? "--Unfiltered_protein_file"
            : ""

    def force_module01_argument =
        force
            ? "--force_genomad --force_checkv"
            : ""

    def interleaved_argument =
        params.interleaved
            .toString()
            .toBoolean()
            ? "--interleaved"
            : ""

    def mapping_delete_argument =
        params.delete_mapping_intermediates
            .toString()
            .toBoolean()
            ? "--delete_files"
            : ""

    def force_mapping_argument =
        force
            ? "--force_read_mapping"
            : ""

    def covered_fraction_argument =
        params.covered_fraction != null
            ? "--covered_fraction ${params.covered_fraction}"
            : ""

    def functional_ads_argument =
        params.functional_ads
            .toString()
            .toBoolean()
            ? "--ADS"
            : ""

    def functional_rdrp_argument =
        params.functional_rdrp
            .toString()
            .toBoolean()
            ? "--RdRP"
            : ""

    def functional_dram_argument =
        params.functional_dram
            .toString()
            .toBoolean()
            ? "--DRAM"
            : ""

    def functional_delete_argument =
        params.delete_functional_intermediates
            .toString()
            .toBoolean()
            ? "--delete_files"
            : ""

    def force_functional_argument =
        force
            ? "--force_prodigal --force_PHROGS " +
              "--force_PFAM --force_outputs"
            : ""

    if (
        force &&
        params.functional_ads
            .toString()
            .toBoolean()
    ) {
        force_functional_argument +=
            " --force_ADS"
    }

    if (
        force &&
        params.functional_rdrp
            .toString()
            .toBoolean()
    ) {
        force_functional_argument +=
            " --force_RdRP"
    }

    def keep_bam_argument =
        params.keep_bam
            .toString()
            .toBoolean()
            ? "--keep_bam"
            : ""

    def binning_delete_argument =
        params.delete_binning_intermediates
            .toString()
            .toBoolean()
            ? "--delete_files"
            : ""

    def binning_sample_group_argument =
        params.binning_sample_group != null
            ? "--binning_sample_group ${params.binning_sample_group}"
            : ""

    def read_mapping_sample_group_argument =
        params.read_mapping_sample_group != null
            ? "--read_mapping_sample_group ${params.read_mapping_sample_group}"
            : ""

    def force_binning_argument =
        force
            ? "--force_vrhyme --force_checkv " +
              "--force_read_mapping --force_outputs"
            : ""

    def summary_force_argument =
        force
            ? "--force"
            : ""

    def miuvig_template_argument =
        params.miuvig_template
            ? "-t ${absPath(params.miuvig_template)}"
            : ""

    """
    set -euo pipefail

    cp "${input_summary}" mvp_input_summary.tsv

    TOOL_ENV="\$(grep '^TOOL_ENV=' "${tools_status}" |
        tail -n 1 |
        cut -d= -f2- || true)"

    if [[ -n "\$TOOL_ENV" && "\$TOOL_ENV" != "SYSTEM" ]]; then
        export PATH="\$TOOL_ENV/bin:\$PATH"
    fi

    if ! command -v mvip >/dev/null 2>&1; then
        echo "ERROR: mvip was not found after tool setup." >&2
        cat "${tools_status}" >&2 || true
        exit 1
    fi

    METADATA="\$(readlink -f "${metadata}")"

    /*
     * Persistent MVP root.
     *
     * Every MVP module receives this exact same directory through -i.
     * MVP therefore creates and reuses its required fixed directories:
     *
     *   00_DATABASES
     *   01_GENOMAD
     *   02_CHECK_V
     *   03_CLUSTERING
     *   04_READ_MAPPING
     *   05_VOTU_TABLES
     *   06_FUNCTIONAL_ANNOTATION
     *   07A_vRHYME_OUTPUT
     *   ...
     *   100_SUMMARIZED_OUTPUTS
     */
    MVP_WORK="${mvp_work}"

    GENOMAD_DB="${genomad_database}"
    CHECKV_DB="${checkv_database}"

    mkdir -p "\$MVP_WORK"
    mkdir -p "${database_root}"
    mkdir -p "\$(dirname "\$GENOMAD_DB")"
    mkdir -p "\$(dirname "\$CHECKV_DB")"

    : > mvp_commands.log
    : > mvp_run.log

    echo "Auxiliary Module 2 MVP started: \$(date)" \
        >> mvp_run.log

    echo "MVP executable: \$(command -v mvip)" \
        >> mvp_run.log

    echo "MVP working directory: \$MVP_WORK" \
        >> mvp_run.log

    echo "MVP metadata: \$METADATA" \
        >> mvp_run.log

    echo "Selected MVP modules: ${ordered_modules.join(',')}" \
        >> mvp_run.log

    echo "geNomad database: \$GENOMAD_DB" \
        >> mvp_run.log

    echo "CheckV database: \$CHECKV_DB" \
        >> mvp_run.log

    echo "Threads: ${task.cpus}" \
        >> mvp_run.log

    echo "----------------------------------------" \
        >> mvp_run.log

    run_command() {
        printf '%q ' "\$@" >> mvp_commands.log
        printf '\\n' >> mvp_commands.log

        echo "==================================================" \
            >> mvp_run.log

        echo "Started: \$(date)" \
            >> mvp_run.log

        printf 'Command: ' >> mvp_run.log
        printf '%q ' "\$@" >> mvp_run.log
        printf '\\n' >> mvp_run.log

        "\$@" >> mvp_run.log 2>&1

        echo "Finished: \$(date)" \
            >> mvp_run.log
    }

    if [[ "${run_module_0}" == "true" ]]; then

        run_command \
            mvip MVP_00_set_up_MVP \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            ${module00_database_flag} \
            ${skip_check_errors_argument} \
            --genomad_db_path "\$GENOMAD_DB" \
            --checkv_db_path "\$CHECKV_DB"
    fi

    if [[ "${run_module_1}" == "true" ]]; then

        run_command \
            mvip MVP_01_run_genomad_checkv \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            --genomad_db_path "\$GENOMAD_DB" \
            --checkv_db_path "\$CHECKV_DB" \
            ${skip_modify_headers_argument} \
            ${min_seq_argument} \
            ${genomad_filter_argument} \
            ${force_module01_argument} \
            --threads ${task.cpus}
    fi

    if [[ "${run_module_2}" == "true" ]]; then

        run_command \
            mvip MVP_02_filter_genomad_checkv \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            --viral_min_genes ${params.viral_min_genes} \
            --host_viral_genes_ratio \
                ${params.host_viral_genes_ratio}
    fi

    if [[ "${run_module_3}" == "true" ]]; then

        run_command \
            mvip MVP_03_do_clustering \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            --min_ani ${params.min_ani} \
            --min_tcov ${params.min_tcov} \
            --min_qcov ${params.min_qcov} \
            --read_type ${params.read_type} \
            ${unfiltered_protein_argument} \
            --threads ${task.cpus}
    fi

    if [[ "${run_module_4}" == "true" ]]; then

        run_command \
            mvip MVP_04_do_read_mapping \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            --read_type ${params.read_type} \
            ${interleaved_argument} \
            ${mapping_delete_argument} \
            ${force_mapping_argument} \
            --threads ${task.cpus}
    fi

    if [[ "${run_module_5}" == "true" ]]; then

        run_command \
            mvip MVP_05_create_vOTU_table \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            ${covered_fraction_argument} \
            --normalization ${params.normalization} \
            --filtration ${params.filtration} \
            --viral_min_genes ${params.viral_min_genes} \
            --host_viral_genes_ratio \
                ${params.host_viral_genes_ratio}
    fi

    if [[ "${run_module_6}" == "true" ]]; then

        run_command \
            mvip MVP_06_do_functional_annotation \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            --fasta_files ${params.functional_fasta_files} \
            --PHROGS_evalue ${params.phrogs_evalue} \
            --PHROGS_score ${params.phrogs_score} \
            --PFAM_evalue ${params.pfam_evalue} \
            --PFAM_score ${params.pfam_score} \
            ${functional_ads_argument} \
            --ADS_evalue ${params.ads_evalue} \
            --ADS_score ${params.ads_score} \
            --ADS_seqid ${params.ads_seqid} \
            ${functional_rdrp_argument} \
            --RdRP_evalue ${params.rdrp_evalue} \
            --RdRP_score ${params.rdrp_score} \
            ${functional_dram_argument} \
            ${functional_delete_argument} \
            ${force_functional_argument} \
            --threads ${task.cpus}
    fi

    if [[ "${run_module_7}" == "true" ]]; then

        run_command \
            mvip MVP_07_do_binning \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            ${binning_sample_group_argument} \
            ${read_mapping_sample_group_argument} \
            ${keep_bam_argument} \
            ${binning_delete_argument} \
            ${interleaved_argument} \
            ${force_binning_argument} \
            --read_type ${params.read_type} \
            --filtration ${params.filtration} \
            --threads ${task.cpus}
    fi

    if [[ "${run_module_99}" == "true" ]]; then

        run_command \
            mvip MVP_99_prep_MIUViG_submission \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            -g "${params.miuvig_identifier}" \
            -s "${params.miuvig_step}" \
            ${miuvig_template_argument}
    fi

    if [[ "${run_module_100}" == "true" ]]; then

        run_command \
            mvip MVP_100_summarize_outputs \
            -i "\$MVP_WORK" \
            -m "\$METADATA" \
            ${summary_force_argument}
    fi

    echo "==================================================" \
        >> mvp_run.log

    echo "Auxiliary Module 2 MVP finished: \$(date)" \
        >> mvp_run.log

    cat > mvp_complete.txt <<EOF
status=complete
finished=\$(date --iso-8601=seconds 2>/dev/null || date)
mvp_work=\$MVP_WORK
metadata=\$METADATA
selected_modules=${ordered_modules.join(',')}
EOF
    """
}