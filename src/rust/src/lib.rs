use std::{path::PathBuf, time::Instant};

#[cfg(feature = "porkchop-integration")]
use std::collections::HashSet;

use bamana::{
    bam::index::IndexKind,
    bam::validate::{
        FindingScope, FindingSeverity, ValidatePayload, ValidationFinding, ValidationMode,
        ValidationSummary,
    },
    commands::check_index::{
        self, CheckIndexPayload, CheckIndexRequest, IndexCompatibility, IndexSupportLevel,
    },
    commands::check_map::{
        self, CheckMapPayload, CheckMapRequest, ConfidenceLevel as MapConfidenceLevel,
        EvidenceSource as MapEvidenceSource, IndexDiagnosticStatus, MappingStatus as MapStatus,
    },
    commands::check_sort::{
        self, CheckSortPayload, CheckSortRequest, ConfidenceLevel as SortConfidenceLevel,
        EvidenceStrength, ObservedOrder,
    },
    commands::check_tag::{
        self, CheckTagMode, CheckTagPayload, CheckTagRequest, CheckTagResult,
        ConfidenceLevel as TagConfidenceLevel,
    },
    commands::summary::{
        self, ConfidenceLevel, FractionSummary, MappingStatus, RecordCountSummary, SummaryMode,
        SummaryPayload, SummaryRequest,
    },
    commands::{
        check_eof::{self, CheckEofRequest, CheckEofResponse},
        validate::{self, ValidateRequest},
        verify::{self, VerifyRequest, VerifyResponse},
    },
    error::AppError,
    formats::probe::{Confidence, ContainerKind, DetectedFormat},
    json::JsonError,
};
use extendr_api::prelude::*;
use pod5_tools::{
    CompareStatus, FilesystemPod5MetadataReader, IntegrityStatus, Pod5CompareReport,
    Pod5FolderInfo, Pod5Manifest, Pod5SubdividePlan, SubdivideStrategy, VerifyCheckStatus,
    VerifyStatus, compare_inputs, find_pod5_directories, folder_info, manifest_from_path,
    read_pod5_file_info, subdivide_plan_from_path, verify_pod5_file,
};
#[cfg(feature = "porkchop-integration")]
use porkchop::{
    cdna::detect_cdna_primer_pair,
    kit::{Kit, KitMetadata, SeqKind, SupportLevel},
    list_supported_kits,
    motif_index::{
        MotifFamily, MotifIndexStrand, cached_motif_index_for_kit, cached_motif_indexes,
    },
};

type SEXP = extendr_api::SEXP;

#[unsafe(no_mangle)]
pub extern "C" fn flounder_rust_capabilities() -> SEXP {
    let porkchop = if cfg!(feature = "porkchop-integration") {
        "|porkchop"
    } else {
        ""
    };
    let payload = format!(
        "flounder.rust_capabilities.v1|flounder-extendr|{}|pod5-tools|bamana{}",
        env!("CARGO_PKG_VERSION"),
        porkchop
    );
    unsafe { Robj::from(payload).get() }
}

#[unsafe(no_mangle)]
pub extern "C" fn flounder_pod5_find(path: SEXP) -> SEXP {
    let result = pod5_find_response(path);
    unsafe { result.get() }
}

#[unsafe(no_mangle)]
pub extern "C" fn flounder_pod5_verify(path: SEXP) -> SEXP {
    let result = pod5_verify_response(path);
    unsafe { result.get() }
}

#[unsafe(no_mangle)]
pub extern "C" fn flounder_pod5_file_info(path: SEXP) -> SEXP {
    let result = pod5_file_info_response(path);
    unsafe { result.get() }
}

#[unsafe(no_mangle)]
pub extern "C" fn flounder_pod5_folder_info(path: SEXP) -> SEXP {
    let result = pod5_folder_info_response(path);
    unsafe { result.get() }
}

#[unsafe(no_mangle)]
pub extern "C" fn flounder_pod5_manifest(path: SEXP) -> SEXP {
    let result = pod5_manifest_response(path);
    unsafe { result.get() }
}

#[unsafe(no_mangle)]
pub extern "C" fn flounder_pod5_compare(left: SEXP, right: SEXP) -> SEXP {
    let result = pod5_compare_response(left, right);
    unsafe { result.get() }
}

#[unsafe(no_mangle)]
pub extern "C" fn flounder_pod5_subdivide_plan(
    path: SEXP,
    strategy: SEXP,
    files_per_chunk: SEXP,
    seconds_per_chunk: SEXP,
    reads_per_chunk: SEXP,
) -> SEXP {
    let result = pod5_subdivide_plan_response(
        path,
        strategy,
        files_per_chunk,
        seconds_per_chunk,
        reads_per_chunk,
    );
    unsafe { result.get() }
}

#[unsafe(no_mangle)]
pub extern "C" fn flounder_bam_summary(
    path: SEXP,
    sample_records: SEXP,
    prefer_index: SEXP,
    include_mapq_hist: SEXP,
    include_flags: SEXP,
    allow_incomplete: SEXP,
) -> SEXP {
    let result = bam_summary_response(
        path,
        sample_records,
        prefer_index,
        include_mapq_hist,
        include_flags,
        allow_incomplete,
    );
    unsafe { result.get() }
}

#[unsafe(no_mangle)]
pub extern "C" fn flounder_bam_verify(path: SEXP) -> SEXP {
    let result = bam_verify_response(path);
    unsafe { result.get() }
}

#[unsafe(no_mangle)]
pub extern "C" fn flounder_bam_validate(
    path: SEXP,
    max_errors: SEXP,
    max_warnings: SEXP,
    header_only: SEXP,
    records: SEXP,
    fail_fast: SEXP,
    include_warnings: SEXP,
) -> SEXP {
    let result = bam_validate_response(
        path,
        max_errors,
        max_warnings,
        header_only,
        records,
        fail_fast,
        include_warnings,
    );
    unsafe { result.get() }
}

#[unsafe(no_mangle)]
pub extern "C" fn flounder_bam_check_eof(path: SEXP) -> SEXP {
    let result = bam_check_eof_response(path);
    unsafe { result.get() }
}

#[unsafe(no_mangle)]
pub extern "C" fn flounder_bam_check_index(path: SEXP, require: SEXP, prefer_csi: SEXP) -> SEXP {
    let result = bam_check_index_response(path, require, prefer_csi);
    unsafe { result.get() }
}

#[unsafe(no_mangle)]
pub extern "C" fn flounder_bam_check_map(
    path: SEXP,
    sample_records: SEXP,
    prefer_index: SEXP,
) -> SEXP {
    let result = bam_check_map_response(path, sample_records, prefer_index);
    unsafe { result.get() }
}

#[unsafe(no_mangle)]
pub extern "C" fn flounder_bam_check_sort(path: SEXP, sample_records: SEXP, strict: SEXP) -> SEXP {
    let result = bam_check_sort_response(path, sample_records, strict);
    unsafe { result.get() }
}

#[unsafe(no_mangle)]
pub extern "C" fn flounder_bam_check_tag(
    path: SEXP,
    tag: SEXP,
    sample_records: SEXP,
    full_scan: SEXP,
    require_type: SEXP,
    count_hits: SEXP,
) -> SEXP {
    let result = bam_check_tag_response(
        path,
        tag,
        sample_records,
        full_scan,
        require_type,
        count_hits,
    );
    unsafe { result.get() }
}

#[unsafe(no_mangle)]
pub extern "C" fn flounder_library_kit_candidates(reads: SEXP, read_ids: SEXP) -> SEXP {
    let result = library_kit_candidates_response(reads, read_ids);
    unsafe { result.get() }
}

#[unsafe(no_mangle)]
pub extern "C" fn flounder_library_adapter_primer_evidence(
    reads: SEXP,
    read_ids: SEXP,
    kit_id: SEXP,
) -> SEXP {
    #[cfg(feature = "porkchop-integration")]
    let result = library_motif_evidence_response(
        reads,
        read_ids,
        kit_id,
        &[MotifFamily::Adapter, MotifFamily::Primer],
    );
    #[cfg(not(feature = "porkchop-integration"))]
    let result = library_motif_evidence_response(reads, read_ids, kit_id, &[]);
    unsafe { result.get() }
}

#[unsafe(no_mangle)]
pub extern "C" fn flounder_library_barcode_evidence(
    reads: SEXP,
    read_ids: SEXP,
    kit_id: SEXP,
) -> SEXP {
    #[cfg(feature = "porkchop-integration")]
    let result = library_motif_evidence_response(
        reads,
        read_ids,
        kit_id,
        &[MotifFamily::Barcode, MotifFamily::Flank],
    );
    #[cfg(not(feature = "porkchop-integration"))]
    let result = library_motif_evidence_response(reads, read_ids, kit_id, &[]);
    unsafe { result.get() }
}

#[unsafe(no_mangle)]
pub extern "C" fn flounder_library_cdna_primer_evidence(
    reads: SEXP,
    read_ids: SEXP,
    kit_id: SEXP,
) -> SEXP {
    let result = library_cdna_primer_evidence_response(reads, read_ids, kit_id);
    unsafe { result.get() }
}

#[unsafe(no_mangle)]
pub extern "C" fn flounder_grammateus_render_figure_html(figure: SEXP) -> SEXP {
    let result = grammateus_render_response(figure, "figure_html");
    unsafe { result.get() }
}

#[unsafe(no_mangle)]
pub extern "C" fn flounder_grammateus_render_figure_pdf(figure: SEXP) -> SEXP {
    let result = grammateus_render_response(figure, "figure_pdf");
    unsafe { result.get() }
}

fn grammateus_render_response(_element: SEXP, operation: &str) -> Robj {
    list!(
        ok = false,
        data = r!(()),
        error = format!(
            "The private Grammateus report renderer is not linked into this floundeR build. \
             Install or build floundeR with an authorized Grammateus runtime before calling {}.",
            operation
        ),
        category = "runtime_unavailable",
        operation = operation
    )
    .into()
}

#[cfg(feature = "porkchop-integration")]
fn library_kit_candidates_response(reads: SEXP, read_ids: SEXP) -> Robj {
    let reads = unsafe { Robj::from_sexp(reads) };
    let read_ids = unsafe { Robj::from_sexp(read_ids) };
    let Ok((reads, read_ids)) = library_read_inputs(&reads, &read_ids) else {
        return library_error_response(
            "argument",
            "`reads` and `read_ids` must be non-missing character vectors of equal length.",
        );
    };

    match library_kit_candidates_data_frame(&reads, &read_ids) {
        Ok(data) => list!(ok = true, data = data, error = r!(()), category = r!(())).into(),
        Err(message) => library_error_response("porkchop", &message),
    }
}

#[cfg(not(feature = "porkchop-integration"))]
fn library_kit_candidates_response(_reads: SEXP, _read_ids: SEXP) -> Robj {
    library_porkchop_unavailable_response()
}

#[cfg(feature = "porkchop-integration")]
fn library_motif_evidence_response(
    reads: SEXP,
    read_ids: SEXP,
    kit_id: SEXP,
    families: &[MotifFamily],
) -> Robj {
    let reads = unsafe { Robj::from_sexp(reads) };
    let read_ids = unsafe { Robj::from_sexp(read_ids) };
    let kit_id = unsafe { Robj::from_sexp(kit_id) };
    let Ok((reads, read_ids)) = library_read_inputs(&reads, &read_ids) else {
        return library_error_response(
            "argument",
            "`reads` and `read_ids` must be non-missing character vectors of equal length.",
        );
    };
    let Some(kit_id) = kit_id.as_str() else {
        return library_error_response(
            "argument",
            "`kit_id` must be a non-missing character scalar.",
        );
    };

    match library_motif_evidence_data_frame(&reads, &read_ids, kit_id, families) {
        Ok(data) => list!(ok = true, data = data, error = r!(()), category = r!(())).into(),
        Err(message) => library_error_response("kit", &message),
    }
}

#[cfg(not(feature = "porkchop-integration"))]
fn library_motif_evidence_response(
    _reads: SEXP,
    _read_ids: SEXP,
    _kit_id: SEXP,
    _families: &[()],
) -> Robj {
    library_porkchop_unavailable_response()
}

#[cfg(feature = "porkchop-integration")]
fn library_cdna_primer_evidence_response(reads: SEXP, read_ids: SEXP, kit_id: SEXP) -> Robj {
    let reads = unsafe { Robj::from_sexp(reads) };
    let read_ids = unsafe { Robj::from_sexp(read_ids) };
    let kit_id = unsafe { Robj::from_sexp(kit_id) };
    let Ok((reads, read_ids)) = library_read_inputs(&reads, &read_ids) else {
        return library_error_response(
            "argument",
            "`reads` and `read_ids` must be non-missing character vectors of equal length.",
        );
    };
    let Some(kit_id) = kit_id.as_str() else {
        return library_error_response(
            "argument",
            "`kit_id` must be a non-missing character scalar.",
        );
    };

    match library_cdna_primer_evidence_data_frame(&reads, &read_ids, kit_id) {
        Ok(data) => list!(ok = true, data = data, error = r!(()), category = r!(())).into(),
        Err(message) => library_error_response("kit", &message),
    }
}

#[cfg(not(feature = "porkchop-integration"))]
fn library_cdna_primer_evidence_response(_reads: SEXP, _read_ids: SEXP, _kit_id: SEXP) -> Robj {
    library_porkchop_unavailable_response()
}

fn pod5_find_response(path: SEXP) -> Robj {
    let path = unsafe { Robj::from_sexp(path) };
    let Some(path) = path.as_str() else {
        return pod5_error_response("`path` must be a non-missing character scalar.");
    };

    match pod5_find_data_frame(path) {
        Ok(data) => list!(ok = true, data = data, error = r!(())).into(),
        Err(error) => pod5_error_response(&error.to_string()),
    }
}

fn pod5_verify_response(path: SEXP) -> Robj {
    let path = unsafe { Robj::from_sexp(path) };
    let Some(path) = path.as_str() else {
        return pod5_error_response_with_category(
            "`path` must be a non-missing character scalar.",
            "path",
        );
    };

    match pod5_verify_data_frame(path) {
        Ok(data) => list!(ok = true, data = data, error = r!(()), category = r!(())).into(),
        Err(error) => pod5_error_response_with_category(&error.to_string(), "path"),
    }
}

fn pod5_file_info_response(path: SEXP) -> Robj {
    let path = unsafe { Robj::from_sexp(path) };
    let Some(path) = path.as_str() else {
        return pod5_error_response_with_category(
            "`path` must be a non-missing character scalar.",
            "path",
        );
    };

    let reader = FilesystemPod5MetadataReader;
    let path = PathBuf::from(path);
    match read_pod5_file_info(&reader, &path) {
        Ok(info) => list!(
            ok = true,
            data = pod5_file_info_data_frame(info),
            error = r!(()),
            category = r!(())
        )
        .into(),
        Err(error) => pod5_error_response_with_category(&error.to_string(), error.category()),
    }
}

fn pod5_folder_info_response(path: SEXP) -> Robj {
    let path = unsafe { Robj::from_sexp(path) };
    let Some(path) = path.as_str() else {
        return pod5_error_response_with_category(
            "`path` must be a non-missing character scalar.",
            "path",
        );
    };

    let reader = FilesystemPod5MetadataReader;
    let path = PathBuf::from(path);
    match folder_info(&path, &reader) {
        Ok(info) => list!(
            ok = true,
            data = pod5_folder_info_data_frame(info),
            error = r!(()),
            category = r!(())
        )
        .into(),
        Err(error) => pod5_error_response_with_category(
            &error.to_string(),
            pod5_tools_error_category(&error.to_string()),
        ),
    }
}

fn pod5_manifest_response(path: SEXP) -> Robj {
    let path = unsafe { Robj::from_sexp(path) };
    let Some(path) = path.as_str() else {
        return pod5_error_response_with_category(
            "`path` must be a non-missing character scalar.",
            "path",
        );
    };

    let path = PathBuf::from(path);
    match manifest_from_path(&path) {
        Ok(manifest) => list!(
            ok = true,
            data = pod5_manifest_data_frame(manifest),
            error = r!(()),
            category = r!(())
        )
        .into(),
        Err(error) => pod5_error_response_with_category(
            &error.to_string(),
            pod5_tools_error_category(&error.to_string()),
        ),
    }
}

fn pod5_compare_response(left: SEXP, right: SEXP) -> Robj {
    let left = unsafe { Robj::from_sexp(left) };
    let right = unsafe { Robj::from_sexp(right) };
    let Some(left) = left.as_str() else {
        return pod5_error_response_with_category(
            "`left` must be a non-missing character scalar.",
            "path",
        );
    };
    let Some(right) = right.as_str() else {
        return pod5_error_response_with_category(
            "`right` must be a non-missing character scalar.",
            "path",
        );
    };

    match compare_inputs(&PathBuf::from(left), &PathBuf::from(right)) {
        Ok(report) => list!(
            ok = true,
            data = pod5_compare_data_frame(report),
            error = r!(()),
            category = r!(())
        )
        .into(),
        Err(error) => pod5_error_response_with_category(
            &error.to_string(),
            pod5_tools_error_category(&error.to_string()),
        ),
    }
}

fn pod5_subdivide_plan_response(
    path: SEXP,
    strategy: SEXP,
    files_per_chunk: SEXP,
    seconds_per_chunk: SEXP,
    reads_per_chunk: SEXP,
) -> Robj {
    let path = unsafe { Robj::from_sexp(path) };
    let strategy = unsafe { Robj::from_sexp(strategy) };
    let files_per_chunk = unsafe { Robj::from_sexp(files_per_chunk) };
    let seconds_per_chunk = unsafe { Robj::from_sexp(seconds_per_chunk) };
    let reads_per_chunk = unsafe { Robj::from_sexp(reads_per_chunk) };

    let Some(path) = path.as_str() else {
        return pod5_error_response_with_category(
            "`path` must be a non-missing character scalar.",
            "path",
        );
    };
    let Some(strategy) = strategy.as_str() else {
        return pod5_error_response_with_category(
            "`strategy` must be a non-missing character scalar.",
            "format",
        );
    };
    let Some(files_per_chunk) = files_per_chunk.as_str() else {
        return pod5_error_response_with_category(
            "`files_per_chunk` must be supplied as text from the R wrapper.",
            "format",
        );
    };
    let Some(seconds_per_chunk) = seconds_per_chunk.as_str() else {
        return pod5_error_response_with_category(
            "`seconds_per_chunk` must be supplied as text from the R wrapper.",
            "format",
        );
    };
    let Some(reads_per_chunk) = reads_per_chunk.as_str() else {
        return pod5_error_response_with_category(
            "`reads_per_chunk` must be supplied as text from the R wrapper.",
            "format",
        );
    };

    let Ok(strategy) = subdivide_strategy_from_label(strategy) else {
        return pod5_error_response_with_category("unsupported subdivision strategy", "format");
    };
    let Ok(files_per_chunk) = files_per_chunk.parse::<u64>() else {
        return pod5_error_response_with_category("files-per-chunk must be an integer", "format");
    };
    let seconds_per_chunk = optional_u64_from_label(seconds_per_chunk);
    let reads_per_chunk = optional_u64_from_label(reads_per_chunk);
    if seconds_per_chunk.is_err() || reads_per_chunk.is_err() {
        return pod5_error_response_with_category(
            "optional subdivision targets must be empty or integer text",
            "format",
        );
    }

    match subdivide_plan_from_path(
        &PathBuf::from(path),
        strategy,
        files_per_chunk,
        seconds_per_chunk.unwrap(),
        reads_per_chunk.unwrap(),
    ) {
        Ok(plan) => list!(
            ok = true,
            data = pod5_subdivide_plan_data_frame(plan),
            error = r!(()),
            category = r!(())
        )
        .into(),
        Err(error) => pod5_error_response_with_category(
            &error.to_string(),
            pod5_tools_error_category(&error.to_string()),
        ),
    }
}

fn bam_summary_response(
    path: SEXP,
    sample_records: SEXP,
    prefer_index: SEXP,
    include_mapq_hist: SEXP,
    include_flags: SEXP,
    allow_incomplete: SEXP,
) -> Robj {
    let path = unsafe { Robj::from_sexp(path) };
    let sample_records = unsafe { Robj::from_sexp(sample_records) };
    let prefer_index = unsafe { Robj::from_sexp(prefer_index) };
    let include_mapq_hist = unsafe { Robj::from_sexp(include_mapq_hist) };
    let include_flags = unsafe { Robj::from_sexp(include_flags) };
    let allow_incomplete = unsafe { Robj::from_sexp(allow_incomplete) };

    let Some(path) = path.as_str() else {
        return bam_error_response(
            "path",
            "`path` must be a non-missing character scalar.",
            None,
            None,
        );
    };
    let Some(sample_records) = sexp_i32(&sample_records) else {
        return bam_error_response(
            "argument",
            "`sample_records` must be supplied as integer text from the R wrapper.",
            None,
            None,
        );
    };
    if sample_records < 0 {
        return bam_error_response(
            "argument",
            "`sample_records` must be zero or a positive integer.",
            None,
            None,
        );
    }
    let Some(prefer_index) = sexp_bool(&prefer_index) else {
        return bam_error_response(
            "argument",
            "`prefer_index` must be TRUE or FALSE.",
            None,
            None,
        );
    };
    let Some(include_mapq_hist) = sexp_bool(&include_mapq_hist) else {
        return bam_error_response(
            "argument",
            "`include_mapq_hist` must be TRUE or FALSE.",
            None,
            None,
        );
    };
    let Some(include_flags) = sexp_bool(&include_flags) else {
        return bam_error_response(
            "argument",
            "`include_flags` must be TRUE or FALSE.",
            None,
            None,
        );
    };
    let Some(allow_incomplete) = sexp_bool(&allow_incomplete) else {
        return bam_error_response(
            "argument",
            "`allow_incomplete` must be TRUE or FALSE.",
            None,
            None,
        );
    };

    let request = SummaryRequest {
        bam: PathBuf::from(path),
        sample_records: sample_records as usize,
        full_scan: sample_records == 0,
        prefer_index,
        include_mapq_hist,
        include_flags,
        regions: Vec::new(),
        live_progress: false,
        allow_incomplete,
    };
    let started = Instant::now();
    let response =
        summary::run(request).with_analysis_wall_seconds(started.elapsed().as_secs_f64());

    if !response.ok {
        let Some(error) = response.error.as_ref() else {
            return bam_error_response(
                "unknown",
                "Bamana summary failed without structured error metadata.",
                None,
                None,
            );
        };
        return bam_error_response(
            &error.code,
            &error.message,
            error.detail.as_deref(),
            error.hint.as_deref(),
        );
    }

    let Some(payload) = response.data.as_ref() else {
        return bam_error_response(
            "unknown",
            "Bamana summary succeeded without a structured payload.",
            None,
            None,
        );
    };

    list!(
        ok = true,
        data = list!(
            status = bam_status_data_frame(&response, payload),
            evidence = bam_evidence_data_frame(payload),
            header = bam_header_data_frame(payload),
            counts = bam_counts_data_frame(payload.counts.as_ref()),
            fractions = bam_fractions_data_frame(payload.fractions.as_ref(), "full_file"),
            fractions_observed =
                bam_fractions_data_frame(payload.fractions_observed.as_ref(), "observed"),
            mapq = bam_mapq_data_frame(payload),
            mapping = bam_mapping_data_frame(payload),
            anomalies = bam_anomalies_data_frame(payload),
            flag_categories = bam_flag_categories_data_frame(payload),
            references = bam_references_data_frame(payload),
            index_derived = bam_index_derived_data_frame(payload),
            mapq_histogram = bam_mapq_histogram_data_frame(payload)
        ),
        error = r!(()),
        category = r!(())
    )
    .into()
}

fn bam_verify_response(path: SEXP) -> Robj {
    let path = unsafe { Robj::from_sexp(path) };
    let Some(path) = path.as_str() else {
        return bam_error_response(
            "path",
            "`path` must be a non-missing character scalar.",
            None,
            None,
        );
    };

    let request = VerifyRequest {
        bam: PathBuf::from(path),
    };
    let started = Instant::now();
    match verify::run(request) {
        Ok(payload) => list!(
            ok = true,
            data =
                bam_verify_data_frame(path, true, started.elapsed().as_secs_f64(), &payload, None),
            error = r!(()),
            category = r!(())
        )
        .into(),
        Err(error) => bam_app_error_response(error),
    }
}

fn bam_validate_response(
    path: SEXP,
    max_errors: SEXP,
    max_warnings: SEXP,
    header_only: SEXP,
    records: SEXP,
    fail_fast: SEXP,
    include_warnings: SEXP,
) -> Robj {
    let path = unsafe { Robj::from_sexp(path) };
    let max_errors = unsafe { Robj::from_sexp(max_errors) };
    let max_warnings = unsafe { Robj::from_sexp(max_warnings) };
    let header_only = unsafe { Robj::from_sexp(header_only) };
    let records = unsafe { Robj::from_sexp(records) };
    let fail_fast = unsafe { Robj::from_sexp(fail_fast) };
    let include_warnings = unsafe { Robj::from_sexp(include_warnings) };

    let Some(path) = path.as_str() else {
        return bam_error_response(
            "path",
            "`path` must be a non-missing character scalar.",
            None,
            None,
        );
    };
    let Some(max_errors) = sexp_i32(&max_errors) else {
        return bam_error_response("argument", "`max_errors` must be an integer.", None, None);
    };
    let Some(max_warnings) = sexp_i32(&max_warnings) else {
        return bam_error_response("argument", "`max_warnings` must be an integer.", None, None);
    };
    let Some(header_only) = sexp_bool(&header_only) else {
        return bam_error_response(
            "argument",
            "`header_only` must be TRUE or FALSE.",
            None,
            None,
        );
    };
    let Some(fail_fast) = sexp_bool(&fail_fast) else {
        return bam_error_response("argument", "`fail_fast` must be TRUE or FALSE.", None, None);
    };
    let Some(include_warnings) = sexp_bool(&include_warnings) else {
        return bam_error_response(
            "argument",
            "`include_warnings` must be TRUE or FALSE.",
            None,
            None,
        );
    };
    if max_errors < 1 || max_warnings < 0 {
        return bam_error_response(
            "argument",
            "`max_errors` must be positive and `max_warnings` must be non-negative.",
            None,
            None,
        );
    }
    let records = match optional_i32_from_label(&records) {
        Ok(records) => records,
        Err(()) => {
            return bam_error_response(
                "argument",
                "`records` must be empty or a positive integer text value.",
                None,
                None,
            );
        }
    };

    let request = ValidateRequest {
        bam: PathBuf::from(path),
        max_errors: max_errors as usize,
        max_warnings: max_warnings as usize,
        header_only,
        records: records.map(|value| value as u64),
        fail_fast,
        include_warnings,
    };
    let started = Instant::now();
    let response =
        validate::run(request).with_analysis_wall_seconds(started.elapsed().as_secs_f64());

    if let Some(payload) = response.data.as_ref() {
        return list!(
            ok = true,
            data = bam_validate_list(&response, payload),
            command_ok = response.ok,
            error = r!(()),
            category = r!(())
        )
        .into();
    }

    let Some(error) = response.error.as_ref() else {
        return bam_error_response(
            "unknown",
            "Bamana validation failed without structured error metadata.",
            None,
            None,
        );
    };
    bam_error_response(
        &error.code,
        &error.message,
        error.detail.as_deref(),
        error.hint.as_deref(),
    )
}

fn bam_check_eof_response(path: SEXP) -> Robj {
    let path = unsafe { Robj::from_sexp(path) };
    let Some(path) = path.as_str() else {
        return bam_error_response(
            "path",
            "`path` must be a non-missing character scalar.",
            None,
            None,
        );
    };

    let request = CheckEofRequest {
        bam: PathBuf::from(path),
    };
    let started = Instant::now();
    match check_eof::run(request) {
        Ok(payload) => list!(
            ok = true,
            data = bam_check_eof_data_frame(
                path,
                true,
                started.elapsed().as_secs_f64(),
                &payload,
                None
            ),
            error = r!(()),
            category = r!(())
        )
        .into(),
        Err(error) => {
            let json_error = error.to_json_error();
            if json_error.code != "truncated_file" {
                return bam_error_response(
                    &json_error.code,
                    &json_error.message,
                    json_error.detail.as_deref(),
                    json_error.hint.as_deref(),
                );
            }
            let data = bam_check_eof_failure_data_frame(
                path,
                started.elapsed().as_secs_f64(),
                &json_error,
            );
            list!(
                ok = true,
                data = data,
                command_ok = false,
                error = r!(()),
                category = r!(())
            )
            .into()
        }
    }
}

fn bam_check_index_response(path: SEXP, require: SEXP, prefer_csi: SEXP) -> Robj {
    let path = unsafe { Robj::from_sexp(path) };
    let require = unsafe { Robj::from_sexp(require) };
    let prefer_csi = unsafe { Robj::from_sexp(prefer_csi) };
    let Some(path) = path.as_str() else {
        return bam_error_response(
            "path",
            "`path` must be a non-missing character scalar.",
            None,
            None,
        );
    };
    let Some(require) = sexp_bool(&require) else {
        return bam_error_response("argument", "`require` must be TRUE or FALSE.", None, None);
    };
    let Some(prefer_csi) = sexp_bool(&prefer_csi) else {
        return bam_error_response(
            "argument",
            "`prefer_csi` must be TRUE or FALSE.",
            None,
            None,
        );
    };

    let request = CheckIndexRequest {
        bam: PathBuf::from(path),
        require,
        prefer_csi,
    };
    let started = Instant::now();
    let response =
        check_index::run(request).with_analysis_wall_seconds(started.elapsed().as_secs_f64());

    if let Some(payload) = response.data.as_ref() {
        return list!(
            ok = true,
            data = bam_check_index_list(&response, payload),
            command_ok = response.ok,
            error = r!(()),
            category = r!(())
        )
        .into();
    }

    let Some(error) = response.error.as_ref() else {
        return bam_error_response(
            "unknown",
            "Bamana index check failed without structured error metadata.",
            None,
            None,
        );
    };
    bam_error_response(
        &error.code,
        &error.message,
        error.detail.as_deref(),
        error.hint.as_deref(),
    )
}

fn bam_check_map_response(path: SEXP, sample_records: SEXP, prefer_index: SEXP) -> Robj {
    let path = unsafe { Robj::from_sexp(path) };
    let sample_records = unsafe { Robj::from_sexp(sample_records) };
    let prefer_index = unsafe { Robj::from_sexp(prefer_index) };
    let Some(path) = path.as_str() else {
        return bam_error_response(
            "path",
            "`path` must be a non-missing character scalar.",
            None,
            None,
        );
    };
    let Some(sample_records) = sexp_i32(&sample_records) else {
        return bam_error_response(
            "argument",
            "`sample_records` must be supplied as an integer.",
            None,
            None,
        );
    };
    if sample_records < 0 {
        return bam_error_response(
            "argument",
            "`sample_records` must be zero or a positive integer.",
            None,
            None,
        );
    }
    let Some(prefer_index) = sexp_bool(&prefer_index) else {
        return bam_error_response(
            "argument",
            "`prefer_index` must be TRUE or FALSE.",
            None,
            None,
        );
    };

    let request = CheckMapRequest {
        bam: PathBuf::from(path),
        sample_records: sample_records as usize,
        full_scan: sample_records == 0,
        prefer_index,
        regions: Vec::new(),
    };
    let started = Instant::now();
    match check_map::run(request) {
        Ok(payload) => list!(
            ok = true,
            data = bam_check_map_list(path, started.elapsed().as_secs_f64(), &payload),
            error = r!(()),
            category = r!(())
        )
        .into(),
        Err(error) => bam_app_error_response(error),
    }
}

fn bam_check_sort_response(path: SEXP, sample_records: SEXP, strict: SEXP) -> Robj {
    let path = unsafe { Robj::from_sexp(path) };
    let sample_records = unsafe { Robj::from_sexp(sample_records) };
    let strict = unsafe { Robj::from_sexp(strict) };
    let Some(path) = path.as_str() else {
        return bam_error_response(
            "path",
            "`path` must be a non-missing character scalar.",
            None,
            None,
        );
    };
    let Some(sample_records) = sexp_i32(&sample_records) else {
        return bam_error_response(
            "argument",
            "`sample_records` must be supplied as an integer.",
            None,
            None,
        );
    };
    if sample_records < 0 {
        return bam_error_response(
            "argument",
            "`sample_records` must be zero or a positive integer.",
            None,
            None,
        );
    }
    let Some(strict) = sexp_bool(&strict) else {
        return bam_error_response("argument", "`strict` must be TRUE or FALSE.", None, None);
    };

    let request = CheckSortRequest {
        bam: PathBuf::from(path),
        sample_records: sample_records as usize,
        strict,
    };
    let started = Instant::now();
    match check_sort::run(request) {
        Ok(payload) => list!(
            ok = true,
            data = bam_check_sort_data_frame(path, started.elapsed().as_secs_f64(), &payload),
            error = r!(()),
            category = r!(())
        )
        .into(),
        Err(error) => bam_app_error_response(error),
    }
}

fn bam_check_tag_response(
    path: SEXP,
    tag: SEXP,
    sample_records: SEXP,
    full_scan: SEXP,
    require_type: SEXP,
    count_hits: SEXP,
) -> Robj {
    let path = unsafe { Robj::from_sexp(path) };
    let tag = unsafe { Robj::from_sexp(tag) };
    let sample_records = unsafe { Robj::from_sexp(sample_records) };
    let full_scan = unsafe { Robj::from_sexp(full_scan) };
    let require_type = unsafe { Robj::from_sexp(require_type) };
    let count_hits = unsafe { Robj::from_sexp(count_hits) };
    let Some(path) = path.as_str() else {
        return bam_error_response(
            "path",
            "`path` must be a non-missing character scalar.",
            None,
            None,
        );
    };
    let Some(tag) = tag.as_str() else {
        return bam_error_response("argument", "`tag` must be a character scalar.", None, None);
    };
    let Some(sample_records) = sexp_i32(&sample_records) else {
        return bam_error_response(
            "argument",
            "`sample_records` must be supplied as an integer.",
            None,
            None,
        );
    };
    if sample_records < 0 {
        return bam_error_response(
            "argument",
            "`sample_records` must be zero or a positive integer.",
            None,
            None,
        );
    }
    let Some(full_scan) = sexp_bool(&full_scan) else {
        return bam_error_response("argument", "`full_scan` must be TRUE or FALSE.", None, None);
    };
    let Some(require_type) = require_type.as_str() else {
        return bam_error_response(
            "argument",
            "`require_type` must be supplied as character text from the R wrapper.",
            None,
            None,
        );
    };
    let Some(count_hits) = sexp_bool(&count_hits) else {
        return bam_error_response(
            "argument",
            "`count_hits` must be TRUE or FALSE.",
            None,
            None,
        );
    };

    let request = CheckTagRequest {
        bam: PathBuf::from(path),
        tag: tag.to_string(),
        sample_records: sample_records as usize,
        full_scan,
        require_type: (!require_type.is_empty()).then(|| require_type.to_string()),
        count_hits,
    };
    let started = Instant::now();
    let response =
        check_tag::run(request).with_analysis_wall_seconds(started.elapsed().as_secs_f64());

    if let Some(payload) = response.data.as_ref() {
        return list!(
            ok = true,
            data = bam_check_tag_data_frame(&response, payload),
            command_ok = response.ok,
            error = r!(()),
            category = r!(())
        )
        .into();
    }

    let Some(error) = response.error.as_ref() else {
        return bam_error_response(
            "unknown",
            "Bamana tag check failed without structured error metadata.",
            None,
            None,
        );
    };
    bam_error_response(
        &error.code,
        &error.message,
        error.detail.as_deref(),
        error.hint.as_deref(),
    )
}

fn pod5_find_data_frame(path: &str) -> Result<Robj, Box<dyn std::error::Error>> {
    let records = find_pod5_directories(PathBuf::from(path))?;

    let paths = records
        .iter()
        .map(|record| record.path.to_string_lossy().into_owned())
        .collect::<Vec<_>>();
    let pod5_file_counts = records
        .iter()
        .map(|record| record.pod5_file_count as f64)
        .collect::<Vec<_>>();
    let total_bytes = records
        .iter()
        .map(|record| record.total_bytes as f64)
        .collect::<Vec<_>>();
    let oldest_modified_utc = records
        .iter()
        .map(|record| record.oldest_modified_utc.clone())
        .collect::<Vec<_>>();
    let newest_modified_utc = records
        .iter()
        .map(|record| record.newest_modified_utc.clone())
        .collect::<Vec<_>>();

    Ok(data_frame!(
        path = paths,
        pod5_file_count = pod5_file_counts,
        total_bytes = total_bytes,
        oldest_modified_utc = oldest_modified_utc,
        newest_modified_utc = newest_modified_utc
    ))
}

fn pod5_verify_data_frame(path: &str) -> Result<Robj, Box<dyn std::error::Error>> {
    let report = verify_pod5_file(&PathBuf::from(path))?;
    let row_count = report.checks.len();
    let paths = vec![report.path.to_string_lossy().into_owned(); row_count];
    let size_bytes = vec![report.size_bytes as f64; row_count];
    let overall_status = vec![verify_status_label(&report.status).to_string(); row_count];
    let checks = report
        .checks
        .iter()
        .map(|check| check.name.clone())
        .collect::<Vec<_>>();
    let categories = report
        .checks
        .iter()
        .map(|check| check.category.clone())
        .collect::<Vec<_>>();
    let statuses = report
        .checks
        .iter()
        .map(|check| verify_check_status_label(&check.status).to_string())
        .collect::<Vec<_>>();
    let details = report
        .checks
        .iter()
        .map(|check| check.detail.clone())
        .collect::<Vec<_>>();

    Ok(data_frame!(
        path = paths,
        size_bytes = size_bytes,
        overall_status = overall_status,
        check = checks,
        category = categories,
        status = statuses,
        detail = details
    ))
}

fn pod5_file_info_data_frame(info: pod5_tools::Pod5FileInfo) -> Robj {
    let (integrity_status, integrity_reason) = integrity_fields(&info.integrity);
    data_frame!(
        path = vec![info.path.to_string_lossy().into_owned()],
        size_bytes = vec![info.size_bytes as f64],
        flow_cell_id = vec![info.flow_cell_id.unwrap_or_default()],
        sequencing_kit = vec![info.sequencing_kit.unwrap_or_default()],
        read_count = vec![
            info.read_count
                .map(|value| value as f64)
                .unwrap_or(f64::NAN)
        ],
        acquisition_start_utc = vec![info.acquisition_start_utc.unwrap_or_default()],
        duration_seconds = vec![info.duration_seconds.unwrap_or(f64::NAN)],
        pod5_version = vec![info.pod5_version.unwrap_or_default()],
        integrity_status = vec![integrity_status.to_string()],
        integrity_reason = vec![integrity_reason]
    )
}

fn pod5_folder_info_data_frame(info: Pod5FolderInfo) -> Robj {
    let (integrity_status, integrity_reason) = integrity_fields(&info.integrity);
    data_frame!(
        path = vec![info.path.to_string_lossy().into_owned()],
        pod5_file_count = vec![info.pod5_file_count as f64],
        total_bytes = vec![info.total_bytes as f64],
        total_reads = vec![
            info.total_reads
                .map(|value| value as f64)
                .unwrap_or(f64::NAN)
        ],
        flow_cell_ids = vec![info.flow_cell_ids.join(",")],
        sequencing_kits = vec![info.sequencing_kits.join(",")],
        acquisition_start_utc = vec![info.acquisition_start_utc.unwrap_or_default()],
        acquisition_end_utc = vec![info.acquisition_end_utc.unwrap_or_default()],
        integrity_status = vec![integrity_status.to_string()],
        integrity_reason = vec![integrity_reason],
        failed_file_count = vec![info.failed_file_count as f64],
        verification_failed_count = vec![info.verification_failed_count as f64],
        duplicate_file_names = vec![info.duplicate_file_names.join(",")],
        warnings = vec![info.warnings.join("; ")]
    )
}

fn pod5_manifest_data_frame(manifest: Pod5Manifest) -> Robj {
    let row_count = manifest.entries.len();
    let schema_versions = vec![manifest.schema_version as f64; row_count];
    let sources = vec![manifest.source.to_string_lossy().into_owned(); row_count];
    let relative_paths = manifest
        .entries
        .iter()
        .map(|entry| entry.relative_path.to_string_lossy().into_owned())
        .collect::<Vec<_>>();
    let paths = manifest
        .entries
        .iter()
        .map(|entry| entry.path.to_string_lossy().into_owned())
        .collect::<Vec<_>>();
    let size_bytes = manifest
        .entries
        .iter()
        .map(|entry| entry.size_bytes as f64)
        .collect::<Vec<_>>();
    let verification_status = manifest
        .entries
        .iter()
        .map(|entry| verify_status_label(&entry.verification_status).to_string())
        .collect::<Vec<_>>();
    let verification_failed_checks = manifest
        .entries
        .iter()
        .map(|entry| entry.verification_failed_checks as f64)
        .collect::<Vec<_>>();

    data_frame!(
        schema_version = schema_versions,
        source = sources,
        relative_path = relative_paths,
        path = paths,
        size_bytes = size_bytes,
        verification_status = verification_status,
        verification_failed_checks = verification_failed_checks
    )
}

fn pod5_compare_data_frame(report: Pod5CompareReport) -> Robj {
    let mut statuses = Vec::<String>::new();
    let mut kinds = Vec::<String>::new();
    let mut relative_paths = Vec::<String>::new();
    let mut left_size_bytes = Vec::<f64>::new();
    let mut right_size_bytes = Vec::<f64>::new();
    let mut left_verification_status = Vec::<String>::new();
    let mut right_verification_status = Vec::<String>::new();
    let status = compare_status_label(&report.status).to_string();

    for path in report.missing_from_right {
        statuses.push(status.clone());
        kinds.push("missing_from_right".to_string());
        relative_paths.push(path.to_string_lossy().into_owned());
        left_size_bytes.push(f64::NAN);
        right_size_bytes.push(f64::NAN);
        left_verification_status.push(String::new());
        right_verification_status.push(String::new());
    }
    for path in report.missing_from_left {
        statuses.push(status.clone());
        kinds.push("missing_from_left".to_string());
        relative_paths.push(path.to_string_lossy().into_owned());
        left_size_bytes.push(f64::NAN);
        right_size_bytes.push(f64::NAN);
        left_verification_status.push(String::new());
        right_verification_status.push(String::new());
    }
    for change in report.changed {
        statuses.push(status.clone());
        kinds.push("changed".to_string());
        relative_paths.push(change.relative_path.to_string_lossy().into_owned());
        left_size_bytes.push(change.left_size_bytes as f64);
        right_size_bytes.push(change.right_size_bytes as f64);
        left_verification_status
            .push(verify_status_label(&change.left_verification_status).to_string());
        right_verification_status
            .push(verify_status_label(&change.right_verification_status).to_string());
    }
    if report.status == CompareStatus::Match {
        statuses.push(status);
        kinds.push("match".to_string());
        relative_paths.push(String::new());
        left_size_bytes.push(f64::NAN);
        right_size_bytes.push(f64::NAN);
        left_verification_status.push(String::new());
        right_verification_status.push(String::new());
    }

    data_frame!(
        status = statuses,
        kind = kinds,
        relative_path = relative_paths,
        left_size_bytes = left_size_bytes,
        right_size_bytes = right_size_bytes,
        left_verification_status = left_verification_status,
        right_verification_status = right_verification_status
    )
}

fn pod5_subdivide_plan_data_frame(plan: Pod5SubdividePlan) -> Robj {
    let warnings = plan.warnings.join("; ");
    let strategy = subdivide_strategy_label(&plan.strategy).to_string();

    if plan.chunks.is_empty() {
        return data_frame!(
            schema_version = vec![plan.schema_version as f64],
            source = vec![plan.source.to_string_lossy().into_owned()],
            strategy = vec![strategy],
            target = vec![plan.target],
            chunk_index = vec![f64::NAN],
            chunk_label = vec![String::new()],
            file_count = vec![0.0],
            total_bytes = vec![0.0],
            read_count = vec![f64::NAN],
            relative_paths = vec![String::new()],
            warnings = vec![warnings]
        );
    }

    let row_count = plan.chunks.len();
    let schema_versions = vec![plan.schema_version as f64; row_count];
    let sources = vec![plan.source.to_string_lossy().into_owned(); row_count];
    let strategies = vec![strategy; row_count];
    let targets = vec![plan.target; row_count];
    let warnings = vec![warnings; row_count];
    let chunk_index = plan
        .chunks
        .iter()
        .map(|chunk| chunk.index as f64)
        .collect::<Vec<_>>();
    let chunk_label = plan
        .chunks
        .iter()
        .map(|chunk| chunk.label.clone())
        .collect::<Vec<_>>();
    let file_count = plan
        .chunks
        .iter()
        .map(|chunk| chunk.file_count as f64)
        .collect::<Vec<_>>();
    let total_bytes = plan
        .chunks
        .iter()
        .map(|chunk| chunk.total_bytes as f64)
        .collect::<Vec<_>>();
    let read_count = plan
        .chunks
        .iter()
        .map(|chunk| {
            chunk
                .read_count
                .map(|value| value as f64)
                .unwrap_or(f64::NAN)
        })
        .collect::<Vec<_>>();
    let relative_paths = plan
        .chunks
        .iter()
        .map(|chunk| {
            chunk
                .relative_paths
                .iter()
                .map(|path| path.to_string_lossy().into_owned())
                .collect::<Vec<_>>()
                .join(",")
        })
        .collect::<Vec<_>>();

    data_frame!(
        schema_version = schema_versions,
        source = sources,
        strategy = strategies,
        target = targets,
        chunk_index = chunk_index,
        chunk_label = chunk_label,
        file_count = file_count,
        total_bytes = total_bytes,
        read_count = read_count,
        relative_paths = relative_paths,
        warnings = warnings
    )
}

fn bam_status_data_frame(
    response: &bamana::json::CommandResponse<SummaryPayload>,
    payload: &SummaryPayload,
) -> Robj {
    data_frame!(
        schema_version = vec![1.0],
        command = vec![response.command.clone()],
        path = vec![response.path.clone().unwrap_or_default()],
        ok = vec![response.ok],
        analysis_wall_seconds = vec![response.analysis_wall_seconds.unwrap_or(f64::NAN)],
        format = vec![payload.format.to_string()],
        mode = vec![summary_mode_label(payload.mode).to_string()],
        confidence = vec![
            payload
                .confidence
                .map(confidence_label)
                .unwrap_or_default()
                .to_string()
        ],
        semantic_note = vec![payload.semantic_note.clone().unwrap_or_default()]
    )
}

fn bam_evidence_data_frame(payload: &SummaryPayload) -> Robj {
    let Some(evidence) = payload.evidence.as_ref() else {
        return data_frame!(
            schema_version = Vec::<f64>::new(),
            header_used = Vec::<bool>::new(),
            index_used = Vec::<bool>::new(),
            records_scanned = Vec::<f64>::new(),
            full_file_scanned = Vec::<bool>::new()
        );
    };

    data_frame!(
        schema_version = vec![1.0],
        header_used = vec![evidence.header_used],
        index_used = vec![evidence.index_used],
        records_scanned = vec![evidence.records_scanned as f64],
        full_file_scanned = vec![evidence.full_file_scanned]
    )
}

fn bam_header_data_frame(payload: &SummaryPayload) -> Robj {
    let Some(header) = payload.header.as_ref() else {
        return data_frame!(
            schema_version = Vec::<f64>::new(),
            references_defined = Vec::<f64>::new(),
            sort_order = Vec::<String>::new(),
            sub_sort_order = Vec::<String>::new(),
            group_order = Vec::<String>::new()
        );
    };

    data_frame!(
        schema_version = vec![1.0],
        references_defined = vec![header.references_defined as f64],
        sort_order = vec![header.sort_order.clone().unwrap_or_default()],
        sub_sort_order = vec![header.sub_sort_order.clone().unwrap_or_default()],
        group_order = vec![header.group_order.clone().unwrap_or_default()]
    )
}

fn bam_counts_data_frame(counts: Option<&RecordCountSummary>) -> Robj {
    let Some(counts) = counts else {
        return data_frame!(
            schema_version = Vec::<f64>::new(),
            records_examined = Vec::<f64>::new(),
            records_total_known = Vec::<f64>::new(),
            mapped_records = Vec::<f64>::new(),
            unmapped_records = Vec::<f64>::new(),
            primary_records = Vec::<f64>::new(),
            secondary_records = Vec::<f64>::new(),
            supplementary_records = Vec::<f64>::new(),
            duplicate_records = Vec::<f64>::new(),
            qc_fail_records = Vec::<f64>::new(),
            paired_records = Vec::<f64>::new(),
            properly_paired_records = Vec::<f64>::new(),
            read1_records = Vec::<f64>::new(),
            read2_records = Vec::<f64>::new()
        );
    };

    data_frame!(
        schema_version = vec![1.0],
        records_examined = vec![counts.records_examined as f64],
        records_total_known = vec![
            counts
                .records_total_known
                .map(|value| value as f64)
                .unwrap_or(f64::NAN)
        ],
        mapped_records = vec![counts.mapped_records as f64],
        unmapped_records = vec![counts.unmapped_records as f64],
        primary_records = vec![counts.primary_records as f64],
        secondary_records = vec![counts.secondary_records as f64],
        supplementary_records = vec![counts.supplementary_records as f64],
        duplicate_records = vec![counts.duplicate_records as f64],
        qc_fail_records = vec![counts.qc_fail_records as f64],
        paired_records = vec![counts.paired_records as f64],
        properly_paired_records = vec![counts.properly_paired_records as f64],
        read1_records = vec![counts.read1_records as f64],
        read2_records = vec![counts.read2_records as f64]
    )
}

fn bam_fractions_data_frame(fractions: Option<&FractionSummary>, scope: &str) -> Robj {
    let Some(fractions) = fractions else {
        return data_frame!(
            schema_version = Vec::<f64>::new(),
            scope = Vec::<String>::new(),
            fraction_mapped = Vec::<f64>::new(),
            fraction_primary = Vec::<f64>::new(),
            fraction_secondary = Vec::<f64>::new(),
            fraction_supplementary = Vec::<f64>::new(),
            fraction_duplicate = Vec::<f64>::new(),
            fraction_qc_fail = Vec::<f64>::new()
        );
    };

    data_frame!(
        schema_version = vec![1.0],
        scope = vec![scope.to_string()],
        fraction_mapped = vec![fractions.fraction_mapped.unwrap_or(f64::NAN)],
        fraction_primary = vec![fractions.fraction_primary.unwrap_or(f64::NAN)],
        fraction_secondary = vec![fractions.fraction_secondary.unwrap_or(f64::NAN)],
        fraction_supplementary = vec![fractions.fraction_supplementary.unwrap_or(f64::NAN)],
        fraction_duplicate = vec![fractions.fraction_duplicate.unwrap_or(f64::NAN)],
        fraction_qc_fail = vec![fractions.fraction_qc_fail.unwrap_or(f64::NAN)]
    )
}

fn bam_mapq_data_frame(payload: &SummaryPayload) -> Robj {
    let Some(mapq) = payload.mapq.as_ref() else {
        return data_frame!(
            schema_version = Vec::<f64>::new(),
            min = Vec::<f64>::new(),
            max = Vec::<f64>::new(),
            mean = Vec::<f64>::new(),
            zero_count = Vec::<f64>::new()
        );
    };

    data_frame!(
        schema_version = vec![1.0],
        min = vec![mapq.min.map(|value| value as f64).unwrap_or(f64::NAN)],
        max = vec![mapq.max.map(|value| value as f64).unwrap_or(f64::NAN)],
        mean = vec![mapq.mean.unwrap_or(f64::NAN)],
        zero_count = vec![mapq.zero_count as f64]
    )
}

fn bam_mapping_data_frame(payload: &SummaryPayload) -> Robj {
    let Some(mapping) = payload.mapping.as_ref() else {
        return data_frame!(
            schema_version = Vec::<f64>::new(),
            status = Vec::<String>::new(),
            references_with_mapped_reads = Vec::<f64>::new(),
            references_with_mapped_reads_observed = Vec::<f64>::new()
        );
    };

    data_frame!(
        schema_version = vec![1.0],
        status = vec![mapping_status_label(mapping.status).to_string()],
        references_with_mapped_reads = vec![
            mapping
                .references_with_mapped_reads
                .map(|value| value as f64)
                .unwrap_or(f64::NAN)
        ],
        references_with_mapped_reads_observed = vec![
            mapping
                .references_with_mapped_reads_observed
                .map(|value| value as f64)
                .unwrap_or(f64::NAN)
        ]
    )
}

fn bam_anomalies_data_frame(payload: &SummaryPayload) -> Robj {
    let Some(anomalies) = payload.anomalies.as_ref() else {
        return data_frame!(
            schema_version = Vec::<f64>::new(),
            contradictory_mapping_state_records = Vec::<f64>::new()
        );
    };

    data_frame!(
        schema_version = vec![1.0],
        contradictory_mapping_state_records =
            vec![anomalies.contradictory_mapping_state_records as f64]
    )
}

fn bam_flag_categories_data_frame(payload: &SummaryPayload) -> Robj {
    let Some(flags) = payload.flag_categories.as_ref() else {
        return data_frame!(
            schema_version = Vec::<f64>::new(),
            paired_records = Vec::<f64>::new(),
            properly_paired_records = Vec::<f64>::new(),
            secondary_records = Vec::<f64>::new(),
            supplementary_records = Vec::<f64>::new(),
            duplicate_records = Vec::<f64>::new(),
            qc_fail_records = Vec::<f64>::new(),
            read1_records = Vec::<f64>::new(),
            read2_records = Vec::<f64>::new(),
            reverse_strand_records = Vec::<f64>::new()
        );
    };

    data_frame!(
        schema_version = vec![1.0],
        paired_records = vec![flags.paired_records as f64],
        properly_paired_records = vec![flags.properly_paired_records as f64],
        secondary_records = vec![flags.secondary_records as f64],
        supplementary_records = vec![flags.supplementary_records as f64],
        duplicate_records = vec![flags.duplicate_records as f64],
        qc_fail_records = vec![flags.qc_fail_records as f64],
        read1_records = vec![flags.read1_records as f64],
        read2_records = vec![flags.read2_records as f64],
        reverse_strand_records = vec![flags.reverse_strand_records as f64]
    )
}

fn bam_references_data_frame(payload: &SummaryPayload) -> Robj {
    let Some(references) = payload.references.as_ref() else {
        return data_frame!(
            schema_version = Vec::<f64>::new(),
            name = Vec::<String>::new(),
            length = Vec::<f64>::new(),
            mapped_reads = Vec::<f64>::new(),
            unmapped_reads = Vec::<f64>::new(),
            observed_mapped = Vec::<String>::new()
        );
    };

    data_frame!(
        schema_version = vec![1.0; references.len()],
        name = references
            .iter()
            .map(|reference| reference.name.clone())
            .collect::<Vec<_>>(),
        length = references
            .iter()
            .map(|reference| reference.length as f64)
            .collect::<Vec<_>>(),
        mapped_reads = references
            .iter()
            .map(|reference| reference
                .mapped_reads
                .map(|value| value as f64)
                .unwrap_or(f64::NAN))
            .collect::<Vec<_>>(),
        unmapped_reads = references
            .iter()
            .map(|reference| reference
                .unmapped_reads
                .map(|value| value as f64)
                .unwrap_or(f64::NAN))
            .collect::<Vec<_>>(),
        observed_mapped = references
            .iter()
            .map(|reference| optional_bool_label(reference.observed_mapped).to_string())
            .collect::<Vec<_>>()
    )
}

fn bam_index_derived_data_frame(payload: &SummaryPayload) -> Robj {
    let Some(index) = payload.index_derived.as_ref() else {
        return data_frame!(
            schema_version = Vec::<f64>::new(),
            present = Vec::<bool>::new(),
            kind = Vec::<String>::new(),
            used = Vec::<bool>::new(),
            total_mapped_reads = Vec::<f64>::new(),
            total_unmapped_reads = Vec::<f64>::new(),
            references_with_mapped_reads = Vec::<f64>::new(),
            note = Vec::<String>::new()
        );
    };

    data_frame!(
        schema_version = vec![1.0],
        present = vec![index.present],
        kind = vec![
            index
                .kind
                .map(index_kind_label)
                .unwrap_or_default()
                .to_string()
        ],
        used = vec![index.used],
        total_mapped_reads = vec![
            index
                .total_mapped_reads
                .map(|value| value as f64)
                .unwrap_or(f64::NAN)
        ],
        total_unmapped_reads = vec![
            index
                .total_unmapped_reads
                .map(|value| value as f64)
                .unwrap_or(f64::NAN)
        ],
        references_with_mapped_reads = vec![
            index
                .references_with_mapped_reads
                .map(|value| value as f64)
                .unwrap_or(f64::NAN)
        ],
        note = vec![index.note.clone().unwrap_or_default()]
    )
}

fn bam_mapq_histogram_data_frame(payload: &SummaryPayload) -> Robj {
    let Some(histogram) = payload
        .mapq
        .as_ref()
        .and_then(|mapq| mapq.histogram.as_ref())
    else {
        return data_frame!(
            schema_version = Vec::<f64>::new(),
            mapq = Vec::<f64>::new(),
            read_count = Vec::<f64>::new()
        );
    };

    data_frame!(
        schema_version = vec![1.0; histogram.len()],
        mapq = histogram
            .keys()
            .map(|value| *value as f64)
            .collect::<Vec<_>>(),
        read_count = histogram
            .values()
            .map(|value| *value as f64)
            .collect::<Vec<_>>()
    )
}

fn bam_verify_data_frame(
    path: &str,
    ok: bool,
    analysis_wall_seconds: f64,
    payload: &VerifyResponse,
    error: Option<&JsonError>,
) -> Robj {
    data_frame!(
        schema_version = vec![1.0],
        command = vec!["verify".to_string()],
        path = vec![path.to_string()],
        ok = vec![ok],
        analysis_wall_seconds = vec![analysis_wall_seconds],
        detected_format = vec![detected_format_label(payload.detected_format).to_string()],
        container = vec![container_kind_label(payload.container).to_string()],
        is_bam = vec![payload.is_bam],
        shallow_verified = vec![payload.shallow_verified],
        deep_validated = vec![payload.deep_validated],
        confidence = vec![probe_confidence_label(payload.confidence).to_string()],
        checks_performed = vec![payload.checks_performed.join(",")],
        semantic_note = vec![payload.semantic_note.clone()],
        error_code = vec![error.map(|error| error.code.clone()).unwrap_or_default()],
        error_message = vec![error.map(|error| error.message.clone()).unwrap_or_default()],
        error_detail = vec![
            error
                .and_then(|error| error.detail.clone())
                .unwrap_or_default()
        ],
        error_hint = vec![
            error
                .and_then(|error| error.hint.clone())
                .unwrap_or_default()
        ]
    )
}

fn bam_validate_list(
    response: &bamana::json::CommandResponse<ValidatePayload>,
    payload: &ValidatePayload,
) -> Robj {
    list!(
        status = bam_validate_status_data_frame(response, payload),
        summary = bam_validate_summary_data_frame(&payload.summary),
        findings = bam_validate_findings_data_frame(&payload.findings),
        error = bam_validate_error_data_frame(response.error.as_ref())
    )
    .into()
}

fn bam_validate_status_data_frame(
    response: &bamana::json::CommandResponse<ValidatePayload>,
    payload: &ValidatePayload,
) -> Robj {
    data_frame!(
        schema_version = vec![1.0],
        command = vec![response.command.clone()],
        path = vec![response.path.clone().unwrap_or_default()],
        ok = vec![response.ok],
        analysis_wall_seconds = vec![response.analysis_wall_seconds.unwrap_or(f64::NAN)],
        format = vec![payload.format.to_string()],
        mode = vec![validation_mode_label(payload.mode).to_string()],
        valid = vec![payload.valid],
        semantic_note = vec![payload.semantic_note.clone()]
    )
}

fn bam_validate_summary_data_frame(summary: &ValidationSummary) -> Robj {
    data_frame!(
        schema_version = vec![1.0],
        header_valid = vec![summary.header_valid],
        records_examined = vec![summary.records_examined as f64],
        full_file_examined = vec![summary.full_file_examined],
        errors = vec![summary.errors as f64],
        warnings = vec![summary.warnings as f64],
        infos = vec![summary.infos as f64]
    )
}

fn bam_validate_findings_data_frame(findings: &[ValidationFinding]) -> Robj {
    data_frame!(
        schema_version = vec![1.0; findings.len()],
        severity = findings
            .iter()
            .map(|finding| finding_severity_label(finding.severity).to_string())
            .collect::<Vec<_>>(),
        scope = findings
            .iter()
            .map(|finding| finding_scope_label(finding.scope).to_string())
            .collect::<Vec<_>>(),
        code = findings
            .iter()
            .map(|finding| finding.code.clone())
            .collect::<Vec<_>>(),
        message = findings
            .iter()
            .map(|finding| finding.message.clone())
            .collect::<Vec<_>>(),
        record_index = findings
            .iter()
            .map(|finding| finding
                .record_index
                .map(|value| value as f64)
                .unwrap_or(f64::NAN))
            .collect::<Vec<_>>(),
        reference_name = findings
            .iter()
            .map(|finding| finding.reference_name.clone().unwrap_or_default())
            .collect::<Vec<_>>(),
        tag = findings
            .iter()
            .map(|finding| finding.tag.clone().unwrap_or_default())
            .collect::<Vec<_>>()
    )
}

fn bam_validate_error_data_frame(error: Option<&JsonError>) -> Robj {
    let Some(error) = error else {
        return data_frame!(
            schema_version = Vec::<f64>::new(),
            code = Vec::<String>::new(),
            message = Vec::<String>::new(),
            detail = Vec::<String>::new(),
            hint = Vec::<String>::new()
        );
    };

    data_frame!(
        schema_version = vec![1.0],
        code = vec![error.code.clone()],
        message = vec![error.message.clone()],
        detail = vec![error.detail.clone().unwrap_or_default()],
        hint = vec![error.hint.clone().unwrap_or_default()]
    )
}

fn bam_check_eof_data_frame(
    path: &str,
    ok: bool,
    analysis_wall_seconds: f64,
    payload: &CheckEofResponse,
    error: Option<&JsonError>,
) -> Robj {
    data_frame!(
        schema_version = vec![1.0],
        command = vec!["check_eof".to_string()],
        path = vec![path.to_string()],
        ok = vec![ok],
        analysis_wall_seconds = vec![analysis_wall_seconds],
        detected_format = vec![detected_format_label(payload.detected_format).to_string()],
        bgzf_eof_present = vec![payload.bgzf_eof_present],
        complete = vec![payload.complete],
        semantic_note = vec![payload.semantic_note.clone()],
        error_code = vec![error.map(|error| error.code.clone()).unwrap_or_default()],
        error_message = vec![error.map(|error| error.message.clone()).unwrap_or_default()],
        error_detail = vec![
            error
                .and_then(|error| error.detail.clone())
                .unwrap_or_default()
        ],
        error_hint = vec![
            error
                .and_then(|error| error.hint.clone())
                .unwrap_or_default()
        ]
    )
}

fn bam_check_eof_failure_data_frame(
    path: &str,
    analysis_wall_seconds: f64,
    error: &JsonError,
) -> Robj {
    data_frame!(
        schema_version = vec![1.0],
        command = vec!["check_eof".to_string()],
        path = vec![path.to_string()],
        ok = vec![false],
        analysis_wall_seconds = vec![analysis_wall_seconds],
        detected_format = vec!["BAM".to_string()],
        bgzf_eof_present = vec![false],
        complete = vec![false],
        semantic_note = vec!["EOF marker evidence was not present; this is a tail-completeness finding and does not by itself describe BAM semantic validity.".to_string()],
        error_code = vec![error.code.clone()],
        error_message = vec![error.message.clone()],
        error_detail = vec![error.detail.clone().unwrap_or_default()],
        error_hint = vec![error.hint.clone().unwrap_or_default()]
    )
}

fn bam_check_index_list(
    response: &bamana::json::CommandResponse<CheckIndexPayload>,
    payload: &CheckIndexPayload,
) -> Robj {
    list!(
        status = bam_check_index_status_data_frame(response, payload),
        index = bam_check_index_inspection_data_frame(payload),
        candidates = bam_check_index_candidates_data_frame(payload),
        error = bam_validate_error_data_frame(response.error.as_ref())
    )
    .into()
}

fn bam_check_index_status_data_frame(
    response: &bamana::json::CommandResponse<CheckIndexPayload>,
    payload: &CheckIndexPayload,
) -> Robj {
    data_frame!(
        schema_version = vec![1.0],
        command = vec![response.command.clone()],
        path = vec![response.path.clone().unwrap_or_default()],
        ok = vec![response.ok],
        analysis_wall_seconds = vec![response.analysis_wall_seconds.unwrap_or(f64::NAN)],
        format = vec![payload.format.to_string()]
    )
}

fn bam_check_index_inspection_data_frame(payload: &CheckIndexPayload) -> Robj {
    let index = &payload.index;
    data_frame!(
        schema_version = vec![1.0],
        present = vec![index.present],
        selected_path = vec![index.selected_path.clone().unwrap_or_default()],
        kind = vec![
            index
                .kind
                .map(index_kind_label)
                .unwrap_or_default()
                .to_string()
        ],
        support_level = vec![index_support_level_label(index.support_level).to_string()],
        usable = vec![index.usable],
        syntactically_valid = vec![optional_bool_label(index.syntactically_valid).to_string()],
        stale = vec![optional_bool_label(index.stale).to_string()],
        bam_newer_than_index = vec![optional_bool_label(index.bam_newer_than_index).to_string()],
        compatibility = vec![index_compatibility_label(index.compatibility).to_string()],
        notes = vec![payload.notes.join("; ")]
    )
}

fn bam_check_index_candidates_data_frame(payload: &CheckIndexPayload) -> Robj {
    data_frame!(
        schema_version = vec![1.0; payload.candidates.len()],
        path = payload
            .candidates
            .iter()
            .map(|candidate| candidate.path.clone())
            .collect::<Vec<_>>(),
        kind = payload
            .candidates
            .iter()
            .map(|candidate| index_kind_label(candidate.kind).to_string())
            .collect::<Vec<_>>(),
        support_level = payload
            .candidates
            .iter()
            .map(|candidate| index_support_level_label(candidate.support_level).to_string())
            .collect::<Vec<_>>(),
        exists = payload
            .candidates
            .iter()
            .map(|candidate| candidate.exists)
            .collect::<Vec<_>>()
    )
}

fn bam_check_map_list(path: &str, analysis_wall_seconds: f64, payload: &CheckMapPayload) -> Robj {
    list!(
        status = bam_check_map_status_data_frame(path, analysis_wall_seconds, payload),
        index = bam_check_map_index_data_frame(payload),
        summary = bam_check_map_summary_data_frame(payload),
        references = bam_check_map_references_data_frame(payload)
    )
    .into()
}

fn bam_check_map_status_data_frame(
    path: &str,
    analysis_wall_seconds: f64,
    payload: &CheckMapPayload,
) -> Robj {
    data_frame!(
        schema_version = vec![1.0],
        command = vec!["check_map".to_string()],
        path = vec![path.to_string()],
        ok = vec![true],
        analysis_wall_seconds = vec![analysis_wall_seconds],
        format = vec![payload.format.to_string()],
        mapping_status = vec![map_check_mapping_status_label(payload.mapping_status).to_string()],
        has_mapped_reads = vec![optional_bool_label(payload.has_mapped_reads).to_string()],
        evidence_source = vec![map_evidence_source_label(payload.evidence_source).to_string()],
        confidence = vec![map_confidence_label(payload.confidence).to_string()],
        semantic_note = vec![payload.semantic_note.clone()]
    )
}

fn bam_check_map_index_data_frame(payload: &CheckMapPayload) -> Robj {
    data_frame!(
        schema_version = vec![1.0],
        present = vec![payload.index.present],
        kind = vec![
            payload
                .index
                .kind
                .map(index_kind_label)
                .unwrap_or_default()
                .to_string()
        ],
        used = vec![payload.index.used],
        diagnostic_status =
            vec![index_diagnostic_status_label(payload.index.diagnostic_status).to_string()],
        diagnostic_detail = vec![payload.index.diagnostic_detail.clone().unwrap_or_default()]
    )
}

fn bam_check_map_summary_data_frame(payload: &CheckMapPayload) -> Robj {
    let summary = &payload.summary;
    data_frame!(
        schema_version = vec![1.0],
        references_defined = vec![summary.references_defined as f64],
        references_with_mapped_reads = vec![
            summary
                .references_with_mapped_reads
                .map(|value| value as f64)
                .unwrap_or(f64::NAN)
        ],
        total_mapped_reads = vec![
            summary
                .total_mapped_reads
                .map(|value| value as f64)
                .unwrap_or(f64::NAN)
        ],
        total_unmapped_reads = vec![
            summary
                .total_unmapped_reads
                .map(|value| value as f64)
                .unwrap_or(f64::NAN)
        ],
        records_examined = vec![
            summary
                .records_examined
                .map(|value| value as f64)
                .unwrap_or(f64::NAN)
        ],
        mapped_records_observed = vec![
            summary
                .mapped_records_observed
                .map(|value| value as f64)
                .unwrap_or(f64::NAN)
        ],
        unmapped_records_observed = vec![
            summary
                .unmapped_records_observed
                .map(|value| value as f64)
                .unwrap_or(f64::NAN)
        ],
        references_with_mapped_reads_observed = vec![
            summary
                .references_with_mapped_reads_observed
                .map(|value| value as f64)
                .unwrap_or(f64::NAN)
        ],
        inconsistent_records_observed = vec![
            summary
                .inconsistent_records_observed
                .map(|value| value as f64)
                .unwrap_or(f64::NAN)
        ],
        region_records_examined = vec![
            summary
                .region_records_examined
                .map(|value| value as f64)
                .unwrap_or(f64::NAN)
        ],
        region_mapped_records_observed = vec![
            summary
                .region_mapped_records_observed
                .map(|value| value as f64)
                .unwrap_or(f64::NAN)
        ],
        region_unmapped_records_observed = vec![
            summary
                .region_unmapped_records_observed
                .map(|value| value as f64)
                .unwrap_or(f64::NAN)
        ]
    )
}

fn bam_check_map_references_data_frame(payload: &CheckMapPayload) -> Robj {
    data_frame!(
        schema_version = vec![1.0; payload.references.len()],
        name = payload
            .references
            .iter()
            .map(|reference| reference.name.clone())
            .collect::<Vec<_>>(),
        length = payload
            .references
            .iter()
            .map(|reference| reference.length as f64)
            .collect::<Vec<_>>(),
        mapped_reads = payload
            .references
            .iter()
            .map(|reference| reference
                .mapped_reads
                .map(|value| value as f64)
                .unwrap_or(f64::NAN))
            .collect::<Vec<_>>(),
        unmapped_reads = payload
            .references
            .iter()
            .map(|reference| reference
                .unmapped_reads
                .map(|value| value as f64)
                .unwrap_or(f64::NAN))
            .collect::<Vec<_>>(),
        observed = payload
            .references
            .iter()
            .map(|reference| optional_bool_label(reference.observed).to_string())
            .collect::<Vec<_>>()
    )
}

fn bam_check_sort_data_frame(
    path: &str,
    analysis_wall_seconds: f64,
    payload: &CheckSortPayload,
) -> Robj {
    let violation = payload.observed_sort.first_violation.as_ref();
    data_frame!(
        schema_version = vec![1.0],
        command = vec!["check_sort".to_string()],
        path = vec![path.to_string()],
        ok = vec![true],
        analysis_wall_seconds = vec![analysis_wall_seconds],
        format = vec![payload.format.to_string()],
        declared_so = vec![payload.declared_sort.so.clone().unwrap_or_default()],
        declared_ss = vec![payload.declared_sort.ss.clone().unwrap_or_default()],
        declared_go = vec![payload.declared_sort.go.clone().unwrap_or_default()],
        observed_order = vec![observed_order_label(payload.observed_sort.order).to_string()],
        observed_sub_order = vec![payload.observed_sort.sub_order.clone().unwrap_or_default()],
        appears_sorted =
            vec![optional_bool_label(payload.observed_sort.appears_sorted).to_string()],
        records_examined = vec![payload.observed_sort.records_examined as f64],
        first_violation_record_index = vec![
            violation
                .map(|violation| violation.record_index as f64)
                .unwrap_or(f64::NAN)
        ],
        first_violation_reason = vec![
            violation
                .map(|violation| violation.reason.clone())
                .unwrap_or_default()
        ],
        evidence_strength =
            vec![evidence_strength_label(payload.observed_sort.evidence_strength).to_string()],
        header_matches_observation =
            vec![optional_bool_label(payload.agreement.header_matches_observation).to_string()],
        confidence = vec![sort_confidence_label(payload.confidence).to_string()],
        semantic_note = vec![payload.semantic_note.clone()]
    )
}

fn bam_check_tag_data_frame(
    response: &bamana::json::CommandResponse<CheckTagPayload>,
    payload: &CheckTagPayload,
) -> Robj {
    data_frame!(
        schema_version = vec![1.0],
        command = vec![response.command.clone()],
        path = vec![response.path.clone().unwrap_or_default()],
        ok = vec![response.ok],
        analysis_wall_seconds = vec![response.analysis_wall_seconds.unwrap_or(f64::NAN)],
        format = vec![payload.format.to_string()],
        tag = vec![payload.tag.clone()],
        required_type = vec![payload.required_type.clone().unwrap_or_default()],
        mode = vec![check_tag_mode_label(payload.mode).to_string()],
        result = vec![check_tag_result_label(payload.result).to_string()],
        tag_found = vec![payload.tag_found],
        records_examined = vec![payload.records_examined as f64],
        records_with_tag = vec![payload.records_with_tag as f64],
        full_file_scanned = vec![payload.full_file_scanned],
        confidence = vec![
            payload
                .confidence
                .map(tag_confidence_label)
                .unwrap_or_default()
                .to_string()
        ],
        semantic_note = vec![payload.semantic_note.clone().unwrap_or_default()],
        error_code = vec![
            response
                .error
                .as_ref()
                .map(|error| error.code.clone())
                .unwrap_or_default()
        ],
        error_message = vec![
            response
                .error
                .as_ref()
                .map(|error| error.message.clone())
                .unwrap_or_default()
        ],
        error_detail = vec![
            response
                .error
                .as_ref()
                .and_then(|error| error.detail.clone())
                .unwrap_or_default()
        ],
        error_hint = vec![
            response
                .error
                .as_ref()
                .and_then(|error| error.hint.clone())
                .unwrap_or_default()
        ]
    )
}

fn summary_mode_label(mode: SummaryMode) -> &'static str {
    match mode {
        SummaryMode::BoundedScan => "bounded_scan",
        SummaryMode::FullScan => "full_scan",
        SummaryMode::Indeterminate => "indeterminate",
    }
}

fn confidence_label(confidence: ConfidenceLevel) -> &'static str {
    match confidence {
        ConfidenceLevel::High => "high",
        ConfidenceLevel::Medium => "medium",
        ConfidenceLevel::Low => "low",
    }
}

fn mapping_status_label(status: MappingStatus) -> &'static str {
    match status {
        MappingStatus::Mapped => "mapped",
        MappingStatus::Unmapped => "unmapped",
        MappingStatus::Indeterminate => "indeterminate",
    }
}

fn index_kind_label(kind: IndexKind) -> &'static str {
    match kind {
        IndexKind::Bai => "BAI",
        IndexKind::Csi => "CSI",
        IndexKind::Gzi => "GZI",
        IndexKind::Unknown => "UNKNOWN",
    }
}

fn index_support_level_label(level: IndexSupportLevel) -> &'static str {
    match level {
        IndexSupportLevel::Absent => "absent",
        IndexSupportLevel::ReadWrite => "read_write",
        IndexSupportLevel::DetectOnly => "detect_only",
        IndexSupportLevel::PlanningSidecar => "planning_sidecar",
        IndexSupportLevel::Unsupported => "unsupported",
    }
}

fn index_compatibility_label(compatibility: IndexCompatibility) -> &'static str {
    match compatibility {
        IndexCompatibility::Plausible => "plausible",
        IndexCompatibility::Absent => "absent",
        IndexCompatibility::Stale => "stale",
        IndexCompatibility::MismatchedOrInvalid => "mismatched_or_invalid",
        IndexCompatibility::DetectedButNotSupported => "detected_but_not_supported",
    }
}

fn map_check_mapping_status_label(status: MapStatus) -> &'static str {
    match status {
        MapStatus::Mapped => "mapped",
        MapStatus::Unmapped => "unmapped",
        MapStatus::Indeterminate => "indeterminate",
    }
}

fn map_evidence_source_label(source: MapEvidenceSource) -> &'static str {
    match source {
        MapEvidenceSource::Index => "index",
        MapEvidenceSource::Scan => "scan",
    }
}

fn index_diagnostic_status_label(status: IndexDiagnosticStatus) -> &'static str {
    match status {
        IndexDiagnosticStatus::NotChecked => "not_checked",
        IndexDiagnosticStatus::Absent => "absent",
        IndexDiagnosticStatus::Usable => "usable",
        IndexDiagnosticStatus::Stale => "stale",
        IndexDiagnosticStatus::Unsupported => "unsupported",
        IndexDiagnosticStatus::Malformed => "malformed",
        IndexDiagnosticStatus::MismatchedReference => "mismatched_reference",
        IndexDiagnosticStatus::Incomplete => "incomplete",
        IndexDiagnosticStatus::Disabled => "disabled",
    }
}

fn map_confidence_label(confidence: MapConfidenceLevel) -> &'static str {
    match confidence {
        MapConfidenceLevel::High => "high",
        MapConfidenceLevel::Medium => "medium",
        MapConfidenceLevel::Low => "low",
    }
}

fn sort_confidence_label(confidence: SortConfidenceLevel) -> &'static str {
    match confidence {
        SortConfidenceLevel::High => "high",
        SortConfidenceLevel::Medium => "medium",
        SortConfidenceLevel::Low => "low",
    }
}

fn observed_order_label(order: ObservedOrder) -> &'static str {
    match order {
        ObservedOrder::Coordinate => "coordinate",
        ObservedOrder::Queryname => "queryname",
        ObservedOrder::Unsorted => "unsorted",
        ObservedOrder::Indeterminate => "indeterminate",
    }
}

fn evidence_strength_label(strength: EvidenceStrength) -> &'static str {
    match strength {
        EvidenceStrength::Strong => "strong",
        EvidenceStrength::Moderate => "moderate",
        EvidenceStrength::Limited => "limited",
    }
}

fn check_tag_mode_label(mode: CheckTagMode) -> &'static str {
    match mode {
        CheckTagMode::BoundedScan => "bounded_scan",
        CheckTagMode::FullScan => "full_scan",
    }
}

fn check_tag_result_label(result: CheckTagResult) -> &'static str {
    match result {
        CheckTagResult::ObservedPresent => "observed_present",
        CheckTagResult::NotFoundInExaminedRecords => "not_found_in_examined_records",
        CheckTagResult::AbsentInFullScan => "absent_in_full_scan",
        CheckTagResult::Indeterminate => "indeterminate",
    }
}

fn tag_confidence_label(confidence: TagConfidenceLevel) -> &'static str {
    match confidence {
        TagConfidenceLevel::High => "high",
        TagConfidenceLevel::Medium => "medium",
        TagConfidenceLevel::Low => "low",
    }
}

fn optional_bool_label(value: Option<bool>) -> &'static str {
    match value {
        Some(true) => "true",
        Some(false) => "false",
        None => "",
    }
}

fn detected_format_label(format: DetectedFormat) -> &'static str {
    match format {
        DetectedFormat::Bam => "BAM",
        DetectedFormat::Sam => "SAM",
        DetectedFormat::Cram => "CRAM",
        DetectedFormat::Fastq => "FASTQ",
        DetectedFormat::FastqGz => "FASTQ.GZ",
        DetectedFormat::Fasta => "FASTA",
        DetectedFormat::Bed => "BED",
        DetectedFormat::Gff => "GFF",
        DetectedFormat::Unknown => "UNKNOWN",
    }
}

fn container_kind_label(container: ContainerKind) -> &'static str {
    match container {
        ContainerKind::Bgzf => "BGZF",
        ContainerKind::Gzip => "GZIP",
        ContainerKind::PlainText => "PLAIN_TEXT",
        ContainerKind::Binary => "BINARY",
        ContainerKind::Unknown => "UNKNOWN",
    }
}

fn probe_confidence_label(confidence: Confidence) -> &'static str {
    match confidence {
        Confidence::High => "high",
        Confidence::Medium => "medium",
        Confidence::Low => "low",
    }
}

fn validation_mode_label(mode: ValidationMode) -> &'static str {
    match mode {
        ValidationMode::HeaderOnly => "header_only",
        ValidationMode::BoundedRecords => "bounded_records",
        ValidationMode::Full => "full",
    }
}

fn finding_severity_label(severity: FindingSeverity) -> &'static str {
    match severity {
        FindingSeverity::Error => "error",
        FindingSeverity::Warning => "warning",
        FindingSeverity::Info => "info",
    }
}

fn finding_scope_label(scope: FindingScope) -> &'static str {
    match scope {
        FindingScope::File => "file",
        FindingScope::Header => "header",
        FindingScope::Record => "record",
        FindingScope::Aux => "aux",
    }
}

fn verify_status_label(status: &VerifyStatus) -> &'static str {
    match status {
        VerifyStatus::Incomplete => "incomplete",
        VerifyStatus::Failed => "failed",
        VerifyStatus::Passed => "passed",
    }
}

fn verify_check_status_label(status: &VerifyCheckStatus) -> &'static str {
    match status {
        VerifyCheckStatus::Passed => "passed",
        VerifyCheckStatus::Failed => "failed",
        VerifyCheckStatus::NotChecked => "not_checked",
    }
}

fn compare_status_label(status: &CompareStatus) -> &'static str {
    match status {
        CompareStatus::Match => "match",
        CompareStatus::Different => "different",
    }
}

fn subdivide_strategy_from_label(label: &str) -> Result<SubdivideStrategy, ()> {
    match label {
        "file-count" => Ok(SubdivideStrategy::FileCount),
        "elapsed-time" => Ok(SubdivideStrategy::ElapsedTime),
        "read-count" => Ok(SubdivideStrategy::ReadCount),
        "sample-label" => Ok(SubdivideStrategy::SampleLabel),
        _ => Err(()),
    }
}

fn subdivide_strategy_label(strategy: &SubdivideStrategy) -> &'static str {
    match strategy {
        SubdivideStrategy::FileCount => "file-count",
        SubdivideStrategy::ElapsedTime => "elapsed-time",
        SubdivideStrategy::ReadCount => "read-count",
        SubdivideStrategy::SampleLabel => "sample-label",
    }
}

fn optional_u64_from_label(label: &str) -> Result<Option<u64>, std::num::ParseIntError> {
    if label.is_empty() {
        Ok(None)
    } else {
        label.parse::<u64>().map(Some)
    }
}

fn optional_i32_from_label(value: &Robj) -> Result<Option<i32>, ()> {
    if value.is_null() {
        return Ok(None);
    }
    let Some(label) = value.as_str() else {
        return Err(());
    };
    if label.is_empty() {
        return Ok(None);
    }
    let parsed = label.parse::<i32>().map_err(|_| ())?;
    if parsed < 1 {
        return Err(());
    }
    Ok(Some(parsed))
}

fn integrity_fields(integrity: &IntegrityStatus) -> (&'static str, String) {
    match integrity {
        IntegrityStatus::Passed => ("passed", String::new()),
        IntegrityStatus::Failed { reason } => ("failed", reason.clone()),
        IntegrityStatus::NotChecked => ("not_checked", String::new()),
        IntegrityStatus::Unavailable { reason } => ("unavailable", reason.clone()),
    }
}

fn pod5_tools_error_category(message: &str) -> &'static str {
    if message.contains("schema") || message.contains("unsupported manifest schema") {
        "schema"
    } else if message.contains("format")
        || message.contains("expects a .pod5 file")
        || message.contains("failed to parse manifest")
    {
        "format"
    } else if message.contains("integrity") {
        "integrity"
    } else if message.contains("failed to inspect")
        || message.contains("failed to read")
        || message.contains("expects a directory")
        || message.contains("expects a file or directory")
    {
        "path"
    } else {
        "unknown"
    }
}

fn sexp_i32(value: &Robj) -> Option<i32> {
    value.as_integer()
}

fn sexp_bool(value: &Robj) -> Option<bool> {
    value.as_bool()
}

fn bam_error_response(code: &str, message: &str, detail: Option<&str>, hint: Option<&str>) -> Robj {
    list!(
        ok = false,
        data = r!(()),
        error = message,
        category = bam_error_category(code),
        code = code,
        detail = detail.unwrap_or_default(),
        hint = hint.unwrap_or_default()
    )
    .into()
}

fn bam_app_error_response(error: AppError) -> Robj {
    let error = error.to_json_error();
    bam_error_response(
        &error.code,
        &error.message,
        error.detail.as_deref(),
        error.hint.as_deref(),
    )
}

#[derive(Debug)]
#[cfg(feature = "porkchop-integration")]
struct KitCandidateRow {
    kit_id: String,
    description: String,
    chemistry: String,
    workflow: String,
    kit_family: String,
    lifecycle_status: String,
    support_level: String,
    legacy: bool,
    introduced_year: f64,
    retired_year: f64,
    score: f64,
    normalized_score: f64,
    matched_motifs: usize,
    total_hits: usize,
    provenance_source: String,
    provenance_appendix: String,
    provenance_notes: String,
    source_urls: String,
    validation_status: String,
    known_limitations: String,
}

#[cfg(feature = "porkchop-integration")]
fn library_read_inputs(reads: &Robj, read_ids: &Robj) -> Result<(Vec<String>, Vec<String>), ()> {
    let Some(reads) = reads.as_str_vector() else {
        return Err(());
    };
    let Some(read_ids) = read_ids.as_str_vector() else {
        return Err(());
    };
    if reads.len() != read_ids.len() {
        return Err(());
    }
    Ok((
        reads.iter().map(|value| (*value).to_string()).collect(),
        read_ids.iter().map(|value| (*value).to_string()).collect(),
    ))
}

#[cfg(feature = "porkchop-integration")]
fn library_kit_candidates_data_frame(
    reads: &[String],
    _read_ids: &[String],
) -> Result<Robj, String> {
    let indexes = cached_motif_indexes().map_err(|error| format!("{error:?}"))?;
    let mut rows = Vec::new();

    for index in indexes {
        let Some(kit) = list_supported_kits()
            .iter()
            .find(|kit| kit.id.0 == index.kit_id().0)
        else {
            continue;
        };
        let mut score = 0.0_f64;
        let mut total_hits = 0_usize;
        let mut matched = HashSet::<String>::new();

        for read in reads {
            for hit in index.find_bidirectional_matches_in(read.as_bytes()) {
                total_hits += 1;
                score += motif_family_weight(hit.entry().family());
                matched.insert(hit.entry().name().to_string());
            }
        }

        let metadata = kit.metadata;
        rows.push(KitCandidateRow {
            kit_id: kit.id.0.to_string(),
            description: kit.description.to_string(),
            chemistry: kit.chemistry.to_string(),
            workflow: metadata.workflow.to_string(),
            kit_family: metadata.family.to_string(),
            lifecycle_status: metadata.status.to_string(),
            support_level: metadata.support_level.to_string(),
            legacy: kit.legacy,
            introduced_year: option_u16_to_f64(metadata.introduced_year),
            retired_year: option_u16_to_f64(metadata.retired_year),
            score,
            normalized_score: 0.0,
            matched_motifs: matched.len(),
            total_hits,
            provenance_source: metadata.provenance.source.to_string(),
            provenance_appendix: metadata.provenance.appendix.unwrap_or_default().to_string(),
            provenance_notes: metadata.provenance.notes.unwrap_or_default().to_string(),
            source_urls: metadata.source_urls.join(";"),
            validation_status: support_level_validation_status(metadata.support_level).to_string(),
            known_limitations: kit_known_limitations(kit, &metadata),
        });
    }

    normalize_candidate_scores(&mut rows);
    rows.sort_by(|left, right| {
        right
            .normalized_score
            .partial_cmp(&left.normalized_score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                right
                    .score
                    .partial_cmp(&left.score)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| left.kit_id.cmp(&right.kit_id))
    });

    Ok(data_frame!(
        schema_version = vec!["flounder.library_kit_candidates.v1".to_string(); rows.len()],
        kit_id = rows
            .iter()
            .map(|row| row.kit_id.clone())
            .collect::<Vec<_>>(),
        description = rows
            .iter()
            .map(|row| row.description.clone())
            .collect::<Vec<_>>(),
        chemistry = rows
            .iter()
            .map(|row| row.chemistry.clone())
            .collect::<Vec<_>>(),
        workflow = rows
            .iter()
            .map(|row| row.workflow.clone())
            .collect::<Vec<_>>(),
        kit_family = rows
            .iter()
            .map(|row| row.kit_family.clone())
            .collect::<Vec<_>>(),
        lifecycle_status = rows
            .iter()
            .map(|row| row.lifecycle_status.clone())
            .collect::<Vec<_>>(),
        support_level = rows
            .iter()
            .map(|row| row.support_level.clone())
            .collect::<Vec<_>>(),
        legacy = rows.iter().map(|row| row.legacy).collect::<Vec<_>>(),
        introduced_year = rows
            .iter()
            .map(|row| row.introduced_year)
            .collect::<Vec<_>>(),
        retired_year = rows.iter().map(|row| row.retired_year).collect::<Vec<_>>(),
        score = rows.iter().map(|row| row.score).collect::<Vec<_>>(),
        normalized_score = rows
            .iter()
            .map(|row| row.normalized_score)
            .collect::<Vec<_>>(),
        score_kind = vec!["heuristic_evidence_score".to_string(); rows.len()],
        matched_motifs = rows
            .iter()
            .map(|row| row.matched_motifs as f64)
            .collect::<Vec<_>>(),
        total_hits = rows
            .iter()
            .map(|row| row.total_hits as f64)
            .collect::<Vec<_>>(),
        provenance_source = rows
            .iter()
            .map(|row| row.provenance_source.clone())
            .collect::<Vec<_>>(),
        provenance_appendix = rows
            .iter()
            .map(|row| row.provenance_appendix.clone())
            .collect::<Vec<_>>(),
        provenance_notes = rows
            .iter()
            .map(|row| row.provenance_notes.clone())
            .collect::<Vec<_>>(),
        source_urls = rows
            .iter()
            .map(|row| row.source_urls.clone())
            .collect::<Vec<_>>(),
        validation_status = rows
            .iter()
            .map(|row| row.validation_status.clone())
            .collect::<Vec<_>>(),
        known_limitations = rows
            .iter()
            .map(|row| row.known_limitations.clone())
            .collect::<Vec<_>>()
    ))
}

#[cfg(feature = "porkchop-integration")]
fn library_motif_evidence_data_frame(
    reads: &[String],
    read_ids: &[String],
    kit_id: &str,
    families: &[MotifFamily],
) -> Result<Robj, String> {
    let Some(index) = cached_motif_index_for_kit(kit_id).map_err(|error| format!("{error:?}"))?
    else {
        return Err(format!("Unknown Porkchop kit id `{kit_id}`."));
    };
    let Some(kit) = list_supported_kits().iter().find(|kit| kit.id.0 == kit_id) else {
        return Err(format!("Unknown Porkchop kit id `{kit_id}`."));
    };

    let mut out_read_ids = Vec::new();
    let mut out_kit_ids = Vec::new();
    let mut motif_names = Vec::new();
    let mut motif_kinds = Vec::new();
    let mut motif_families = Vec::new();
    let mut motif_sources = Vec::new();
    let mut strands = Vec::new();
    let mut starts = Vec::new();
    let mut ends = Vec::new();
    let mut match_semantics = Vec::new();
    let mut edit_distances = Vec::new();
    let mut provenance_sources = Vec::new();
    let mut support_levels = Vec::new();

    for (read, read_id) in reads.iter().zip(read_ids.iter()) {
        for hit in index.find_bidirectional_matches_in(read.as_bytes()) {
            let entry = hit.entry();
            if !families.contains(&entry.family()) {
                continue;
            }
            let record = entry.record();
            out_read_ids.push(read_id.clone());
            out_kit_ids.push(kit_id.to_string());
            motif_names.push(entry.name().to_string());
            motif_kinds.push(seq_kind_label(entry.kind()).to_string());
            motif_families.push(motif_family_label(entry.family()).to_string());
            motif_sources.push(motif_source_label(entry.source()).to_string());
            strands.push(motif_strand_label(hit.strand()).to_string());
            starts.push(hit.start() as f64);
            ends.push(hit.end() as f64);
            match_semantics.push("ambiguity_aware_exact".to_string());
            edit_distances.push(f64::NAN);
            provenance_sources.push(record.provenance.source.to_string());
            support_levels.push(kit.metadata.support_level.to_string());
        }
    }

    Ok(data_frame!(
        schema_version = vec!["flounder.library_motif_evidence.v1".to_string(); out_read_ids.len()],
        read_id = out_read_ids,
        kit_id = out_kit_ids,
        motif_name = motif_names,
        motif_kind = motif_kinds,
        motif_family = motif_families,
        motif_source = motif_sources,
        strand = strands,
        start = starts,
        end = ends,
        match_semantics = match_semantics,
        edit_distance = edit_distances,
        provenance_source = provenance_sources,
        support_level = support_levels
    ))
}

#[cfg(feature = "porkchop-integration")]
fn library_cdna_primer_evidence_data_frame(
    reads: &[String],
    read_ids: &[String],
    kit_id: &str,
) -> Result<Robj, String> {
    if porkchop::cdna::orientation_rules_for_kit(kit_id).is_none() {
        return Err(format!(
            "Porkchop kit `{kit_id}` has no cDNA orientation rules."
        ));
    }

    let mut out_read_ids = Vec::new();
    let mut kit_ids = Vec::new();
    let mut classes = Vec::new();
    let mut classified = Vec::new();
    let mut full_length = Vec::new();
    let mut five_prime_names = Vec::new();
    let mut five_prime_starts = Vec::new();
    let mut five_prime_ends = Vec::new();
    let mut three_prime_names = Vec::new();
    let mut three_prime_starts = Vec::new();
    let mut three_prime_ends = Vec::new();
    let mut primer_hit_counts = Vec::new();
    let mut workflow_support_notes = Vec::new();
    let mut known_limitations = Vec::new();

    let kit = list_supported_kits().iter().find(|kit| kit.id.0 == kit_id);
    let limitations = kit
        .map(|kit| kit_known_limitations(kit, &kit.metadata))
        .unwrap_or_default();
    let support_note = "primer-pair classification only; rescue, UMI-aware barcode handling, and public-dataset validation are separate claims";

    for (read, read_id) in reads.iter().zip(read_ids.iter()) {
        let Some(result) = detect_cdna_primer_pair(kit_id, read.as_bytes()) else {
            return Err(format!(
                "Porkchop kit `{kit_id}` has no cDNA registry sequence records."
            ));
        };
        out_read_ids.push(read_id.clone());
        kit_ids.push(result.kit_id.clone());
        classes.push(result.class.as_str().to_string());
        classified.push(result.class.is_classified());
        full_length.push(result.class.is_full_length());
        five_prime_names.push(
            result
                .five_prime
                .as_ref()
                .map(|hit| hit.name.clone())
                .unwrap_or_default(),
        );
        five_prime_starts.push(
            result
                .five_prime
                .as_ref()
                .map(|hit| hit.start as f64)
                .unwrap_or(f64::NAN),
        );
        five_prime_ends.push(
            result
                .five_prime
                .as_ref()
                .map(|hit| hit.end as f64)
                .unwrap_or(f64::NAN),
        );
        three_prime_names.push(
            result
                .three_prime
                .as_ref()
                .map(|hit| hit.name.clone())
                .unwrap_or_default(),
        );
        three_prime_starts.push(
            result
                .three_prime
                .as_ref()
                .map(|hit| hit.start as f64)
                .unwrap_or(f64::NAN),
        );
        three_prime_ends.push(
            result
                .three_prime
                .as_ref()
                .map(|hit| hit.end as f64)
                .unwrap_or(f64::NAN),
        );
        primer_hit_counts.push(result.primer_hits.len() as f64);
        workflow_support_notes.push(support_note.to_string());
        known_limitations.push(limitations.clone());
    }

    Ok(data_frame!(
        schema_version =
            vec!["flounder.library_cdna_primer_evidence.v1".to_string(); out_read_ids.len()],
        kit_id = kit_ids,
        read_id = out_read_ids,
        class = classes,
        classified = classified,
        full_length = full_length,
        five_prime_name = five_prime_names,
        five_prime_start = five_prime_starts,
        five_prime_end = five_prime_ends,
        three_prime_name = three_prime_names,
        three_prime_start = three_prime_starts,
        three_prime_end = three_prime_ends,
        primer_hit_count = primer_hit_counts,
        workflow_support_note = workflow_support_notes,
        known_limitations = known_limitations
    ))
}

#[cfg(feature = "porkchop-integration")]
fn normalize_candidate_scores(rows: &mut [KitCandidateRow]) {
    if rows.is_empty() {
        return;
    }
    let max_score = rows
        .iter()
        .map(|row| row.score)
        .fold(f64::NEG_INFINITY, f64::max);
    let exps = rows
        .iter()
        .map(|row| (row.score - max_score).exp())
        .collect::<Vec<_>>();
    let denominator = exps.iter().sum::<f64>().max(1e-12);
    for (row, exp) in rows.iter_mut().zip(exps.iter()) {
        row.normalized_score = exp / denominator;
    }
}

#[cfg(feature = "porkchop-integration")]
fn motif_family_weight(family: MotifFamily) -> f64 {
    match family {
        MotifFamily::Adapter => 3.0,
        MotifFamily::Primer => 2.0,
        MotifFamily::Barcode => 1.0,
        MotifFamily::Flank => 0.5,
    }
}

#[cfg(feature = "porkchop-integration")]
fn option_u16_to_f64(value: Option<u16>) -> f64 {
    value.map(|value| value as f64).unwrap_or(f64::NAN)
}

#[cfg(feature = "porkchop-integration")]
fn support_level_validation_status(level: SupportLevel) -> &'static str {
    match level {
        SupportLevel::Full => "registry_supported_validation_separate",
        SupportLevel::BestEffort => "best_effort_registry_support",
        SupportLevel::Partial => "partial_registry_support",
        SupportLevel::Experimental => "experimental_registry_support",
    }
}

#[cfg(feature = "porkchop-integration")]
fn kit_known_limitations(kit: &Kit, metadata: &KitMetadata) -> String {
    let mut notes = Vec::new();
    if kit.legacy {
        notes.push("legacy chemistry support is conservative");
    }
    match metadata.support_level {
        SupportLevel::Full => {
            notes.push("workflow validation remains separate from registry support")
        }
        SupportLevel::BestEffort => {
            notes.push("best-effort support must not be presented as validated workflow support")
        }
        SupportLevel::Partial => notes
            .push("partial kit model; do not treat as complete trimming or demultiplexing support"),
        SupportLevel::Experimental => {
            notes.push("experimental support; review before production use")
        }
    }
    if matches!(
        metadata.workflow,
        porkchop::kit::Workflow::PCRcDNASequencing | porkchop::kit::Workflow::PCRcDNABarcoding
    ) {
        notes.push("cDNA rescue, UMI-aware barcode handling, and public-dataset validation are separate claims");
    }
    notes.join("; ")
}

#[cfg(feature = "porkchop-integration")]
fn seq_kind_label(kind: SeqKind) -> &'static str {
    match kind {
        SeqKind::AdapterTop => "adapter_top",
        SeqKind::AdapterBottom => "adapter_bottom",
        SeqKind::Primer => "primer",
        SeqKind::Barcode => "barcode",
        SeqKind::Flank => "flank",
    }
}

#[cfg(feature = "porkchop-integration")]
fn motif_family_label(family: MotifFamily) -> &'static str {
    match family {
        MotifFamily::Adapter => "adapter",
        MotifFamily::Primer => "primer",
        MotifFamily::Barcode => "barcode",
        MotifFamily::Flank => "flank",
    }
}

#[cfg(feature = "porkchop-integration")]
fn motif_source_label(source: porkchop::motif_index::MotifSource) -> &'static str {
    match source {
        porkchop::motif_index::MotifSource::AdapterOrPrimer => "adapter_or_primer",
        porkchop::motif_index::MotifSource::BarcodeOrFlank => "barcode_or_flank",
    }
}

#[cfg(feature = "porkchop-integration")]
fn motif_strand_label(strand: MotifIndexStrand) -> &'static str {
    match strand {
        MotifIndexStrand::Forward => "forward",
        MotifIndexStrand::Reverse => "reverse",
    }
}

#[cfg(not(feature = "porkchop-integration"))]
fn library_porkchop_unavailable_response() -> Robj {
    library_error_response(
        "porkchop_unavailable",
        "Porkchop-backed library-preparation evidence is not compiled into this floundeR build. Reinstall from source with CARGO_FEATURE_ARGS=--features=porkchop-integration and a local ../porkchop checkout, or use a build that links the public Porkchop crate.",
    )
}

fn library_error_response(category: &str, message: &str) -> Robj {
    list!(
        ok = false,
        data = r!(()),
        error = message,
        category = category
    )
    .into()
}

fn bam_error_category(code: &str) -> &'static str {
    match code {
        "file_not_found" | "permission_denied" | "io_error" | "unknown_format" => "path",
        "invalid_bam"
        | "invalid_header"
        | "invalid_record"
        | "not_bam"
        | "truncated_file"
        | "unsupported_input_format" => "format",
        "invalid_index" | "missing_index" | "unsupported_index" => "index",
        "argument" | "invalid_tag" | "invalid_tag_type" => "argument",
        _ => "unknown",
    }
}

fn pod5_error_response(message: &str) -> Robj {
    pod5_error_response_with_category(message, "unknown")
}

fn pod5_error_response_with_category(message: &str, category: &str) -> Robj {
    list!(
        ok = false,
        data = r!(()),
        error = message,
        category = category
    )
    .into()
}
