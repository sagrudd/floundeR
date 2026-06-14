use std::path::PathBuf;

use extendr_api::prelude::*;
use pod5_tools::{
    CompareStatus, FilesystemPod5MetadataReader, IntegrityStatus, Pod5CompareReport,
    Pod5FolderInfo, Pod5Manifest, VerifyCheckStatus, VerifyStatus, compare_inputs,
    find_pod5_directories, folder_info, manifest_from_path, read_pod5_file_info, verify_pod5_file,
};

type SEXP = extendr_api::SEXP;

#[unsafe(no_mangle)]
pub extern "C" fn flounder_rust_capabilities() -> SEXP {
    let payload = format!(
        "flounder.rust_capabilities.v1|flounder-extendr|{}|pod5-tools",
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
