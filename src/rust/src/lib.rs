#![no_std]

use core::ffi::{c_char, c_void};
use core::panic::PanicInfo;

type SEXP = *mut c_void;

unsafe extern "C" {
    fn Rf_mkString(value: *const c_char) -> SEXP;
}

#[unsafe(no_mangle)]
pub extern "C" fn flounder_rust_capabilities() -> SEXP {
    unsafe {
        const PAYLOAD: &[u8] = concat!(
            "flounder.rust_capabilities.v1|flounder-extendr|",
            env!("CARGO_PKG_VERSION"),
            "\0"
        )
        .as_bytes();
        Rf_mkString(PAYLOAD.as_ptr().cast())
    }
}

#[panic_handler]
fn panic(_info: &PanicInfo) -> ! {
    loop {}
}
