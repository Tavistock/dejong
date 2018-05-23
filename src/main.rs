extern crate image;

use std::fs::File;
use std::usize;
use std::f64;

struct DeJong (f64, f64, f64, f64);

// fn scale_linear(v: f64) {255.0 * v}

fn scale(v: f64) -> f64 {v.sqrt()}
// fn scale(v: f64, t: u8) -> f64 {(t as f64 / 2f64.ln()) * (v + 1f64).ln()}
// fn scale(v: f64, s: f64, t: u8) -> f64 {
//     (t as f64 / (2f64 + s).ln()) * ((v * s) + 1f64).ln()
// }

fn lerp(bg: &image::Rgb<u8>, fg: &image::Rgb<u8>, t: f64) -> image::Rgb<u8> {
    let max_t = u8::max_value() as f64;
    let (bg_r, bg_g, bg_b) = (bg.data[0], bg.data[1], bg.data[2]);
    let (fg_r, fg_g, fg_b) = (fg.data[0], fg.data[1], fg.data[2]);
    let (bg_r, bg_g, bg_b) = (bg_r as f64 / max_t,
                              bg_g as f64 / max_t,
                              bg_b as f64 / max_t);
    let (fg_r, fg_g, fg_b) = (fg_r as f64 / max_t,
                              fg_g as f64 / max_t,
                              fg_b as f64 / max_t);
    let (out_r, out_g, out_b) = (bg_r + t * (fg_r - bg_r),
                                 bg_g + t * (fg_g - bg_g),
                                 bg_b + t * (fg_b - bg_b));
    image::Rgb([
        (max_t * out_r) as u8,
        (max_t * out_g) as u8,
        (max_t * out_b) as u8,
    ])
}

// fn next_point(x: f64, y: f64, params: &DeJong) -> (f64, f64){
//     let &DeJong (a, b, c, d) = params;
//     let nx = (a * y).sin() - (b * x).cos();
//     let ny = (c * x).sin() - (d * y).cos();
//     (nx, ny)
// }

fn next_point(x: f64, y: f64, params: &DeJong) -> (f64, f64){
    let &DeJong (a, b, c, d) = params;
    let nx = (a * y).sin() + c * (a * x).cos();
    let ny = (b * x).sin() - d * (b * y).cos();
    (nx, ny)
}


// fn next_point(x: f64, y: f64, params: &DeJong) -> (f64, f64){
//     let &DeJong (a, b, c, d) = params;
//     let nx = d * (a * x).sin() - (b *y).sin();
//     let ny = c * (a * x).cos() + (b * y).cos();
//     (nx, ny)
// }

fn clamp(x: f64, min: f64, max: f64) -> f64 {
    if x < min {min} else if x > max {max} else {x}
}

fn intensity(dx: f64, dy: f64) -> f64 {
    clamp((dx.abs() + dy.abs()) / 4., 0., 1.)
}

fn expose(sx: f64, sy: f64) -> (f64, f64, f64, f64) {
    // the floor of the values, to use as the array indexes
    let bx = sx as u32;
    let by = sy as u32;

    // The remainder
    let xr = sx - bx as f64;
    let yr = sy - by as f64;

    // 1 - remainder
    let xr1 = 1.0 - xr;
    let yr1 = 1.0 - yr;

    // calculate the partial pixel intensities
    (1. - (xr  * xr  + yr  * yr ) / 2.,
     1. - (xr  * xr  + yr1 * yr1) / 2.,
     1. - (xr1 * xr1 + yr  * yr ) / 2.,
     1. - (xr1 * xr1 + yr1 * yr1) / 2.)
image::Rgb<u8>}

fn dejong(dim: u32, params: DeJong, len: u32) -> Vec<f64> {
    let mut arr = vec![0f64; (dim * dim) as usize];
    let mut x = 0f64;
    let mut y = 0f64;
    let d2 = dim as f64 / 2.;
    let d5 = dim as f64 / 4. - 1.;
    for _ in 0..len {
        let (nx, ny) = next_point(x, y, &params);
        // scale to fit the dimennsions of the image
        let i = intensity(x - nx, y - ny);
        let sx = nx * d5 + d2;
        let sy = ny * d5 + d2;
        let (e, f, g, h) = expose(sx, sy);
        let bx = sx as u32;
        let by = sy as u32;
        x = nx;
        y = ny;
        arr[((by  ) * dim + bx  ) as usize] += e * i;
        arr[((by  ) * dim + bx+1) as usize] += f * i;
        arr[((by+1) * dim + bx  ) as usize] += g * i;
        arr[((by+1) * dim + bx+1) as usize] += h * i;
    }
    arr
}

fn scale_dejong(mut vec: Vec<f64>) -> Vec<f64> {
    let max = vec.iter().cloned().fold(f64::NAN, f64::max);
    vec.iter_mut()
        .map(|n| scale(*n * (1./ max)))
        .collect::<Vec<_>>()
}

fn color_dejong(mut vec: Vec<f64>, c2: image::Rgb<u8>, c1: image::Rgb<u8>) ->
    Vec<image::Rgb<u8>>
{
    vec.iter_mut()
        .map(|n| lerp(&c1, &c2, *n))
        .collect::<Vec<image::Rgb<u8>>>()
}

fn main() {
    let iterations = 5_000_000;
    // let params = DeJong(-2.24, 0.43, -0.65, -2.43);
    // let params = DeJong( 2.01, -2.53, 1.61, -0.33);
    // let params = DeJong(-2., -2., -1.2, 2.);
    // let params = DeJong(1.4, -2.3, 2.4, -2.1);
    // let params = DeJong(-2., 0.43, -1.65, 2.);
    // let params = DeJong(1.40, 1.56, 1.40, -1.);
    let params = DeJong(1., -1.8, 1.6, 2.);
    let dim = 800;
    let c1 = image::Rgb([247, 202, 201]);
    let c2 = image::Rgb([23, 23, 23]);
    let mut imgbuf = image::ImageBuffer::new(dim, dim);
    let arr = color_dejong(scale_dejong(dejong(dim, params, iterations)), c1, c2);
    for (i, pixel) in arr.iter().enumerate() {
        let y = i as u32 % dim;
        let x = i as u32 / dim;
        *imgbuf.get_pixel_mut(x, y) = *pixel;
    }
    let ref mut fout = File::create("fractal.png").unwrap();
    image::ImageRgb8(imgbuf).save(fout, image::PNG).unwrap();
}
