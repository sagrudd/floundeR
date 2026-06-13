#![no_std]

use core::ffi::{c_char, c_void};
use core::panic::PanicInfo;

type SEXP = *mut c_void;

extern "C" {
    fn Rf_mkString(value: *const c_char) -> SEXP;
}

#[no_mangle]
pub extern "C" fn flounder_rust_capabilities() -> SEXP {
    unsafe {
        const PAYLOAD: &[u8] = b"flounder.rust_capabilities.v1|flounder-extendr|0.1.6\0";
        Rf_mkString(PAYLOAD.as_ptr().cast())
    }
}

#[panic_handler]
fn panic(_info: &PanicInfo) -> ! {
    loop {}
}
