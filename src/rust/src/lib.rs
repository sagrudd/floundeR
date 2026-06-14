use std::{path::PathBuf, time::Instant};

use bamana::{
    bam::index::IndexKind,
    commands::summary::{
        self, ConfidenceLevel, FractionSummary, MappingStatus, RecordCountSummary, SummaryMode,
        SummaryPayload, SummaryRequest,
    },
};
use extendr_api::prelude::*;
use pod5_tools::{
    CompareStatus, FilesystemPod5MetadataReader, IntegrityStatus, Pod5CompareReport,
    Pod5FolderInfo, Pod5Manifest, Pod5SubdividePlan, SubdivideStrategy, VerifyCheckStatus,
    VerifyStatus, compare_inputs, find_pod5_directories, folder_info, manifest_from_path,
    read_pod5_file_info, subdivide_plan_from_path, verify_pod5_file,
};

type SEXP = extendr_api::SEXP;

#[unsafe(no_mangle)]
pub extern "C" fn flounder_rust_capabilities() -> SEXP {
    let payload = format!(
        "flounder.rust_capabilities.v1|flounder-extendr|{}|pod5-tools|bamana",
        env!("CARGO_PKG_VERSION")
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

fn optional_bool_label(value: Option<bool>) -> &'static str {
    match value {
        Some(true) => "true",
        Some(false) => "false",
        None => "",
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

fn bam_error_category(code: &str) -> &'static str {
    match code {
        "file_not_found" | "permission_denied" | "io_error" | "unknown_format" => "path",
        "invalid_bam"
        | "invalid_header"
        | "invalid_record"
        | "not_bam"
        | "unsupported_input_format" => "format",
        "invalid_index" | "missing_index" | "unsupported_index" => "index",
        "argument" => "argument",
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
