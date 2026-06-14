use std::path::PathBuf;

use extendr_api::prelude::*;
use pod5_tools::{
    FilesystemPod5MetadataReader, IntegrityStatus, VerifyCheckStatus, VerifyStatus,
    find_pod5_directories, read_pod5_file_info, verify_pod5_file,
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

fn integrity_fields(integrity: &IntegrityStatus) -> (&'static str, String) {
    match integrity {
        IntegrityStatus::Passed => ("passed", String::new()),
        IntegrityStatus::Failed { reason } => ("failed", reason.clone()),
        IntegrityStatus::NotChecked => ("not_checked", String::new()),
        IntegrityStatus::Unavailable { reason } => ("unavailable", reason.clone()),
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
