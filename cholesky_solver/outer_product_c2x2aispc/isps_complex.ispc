struct DComplex {
  double re, im;
};
inline
uniform DComplex conj(uniform  DComplex a) {
    uniform DComplex result;
    result.re =  a.re;
    result.im = -a.im;
    return result;
}
inline
uniform double norm(uniform  DComplex a) {
    return a.re*a.re + a.im*a.im;
}
inline
uniform DComplex operator*(uniform  DComplex a, uniform DComplex b) {
    uniform DComplex result;
    result.re = a.re * b.re - a.im * b.im;
    result.im = a.re * b.im + a.im * b.re;
    return result;
}
inline
uniform DComplex operator-(uniform DComplex a, uniform DComplex b) {
    uniform DComplex result;
    result.re = a.re - b.re;
    result.im = a.im - b.im;
    return result;
}
inline
uniform struct DComplex operator*(uniform struct DComplex a, uniform double b) {
    uniform struct DComplex result;
    result.re = a.re * b;
    result.im = a.im * b;
    return result;
}
