use polynomial::*;

fn main() {
    let point_ys = vec![Fq::from(0), Fq::from(4), Fq::from(16)];
    let point_xs = vec![Fq::from(0), Fq::from(2), Fq::from(4)];

    let poly = DenseUnivariatePolynomial::interpolate(point_ys, point_xs);
    // dbg!(&poly);
    // let poly = DenseUnivariatePolynomial::interpolate(point_ys, point_xs);
    // assert_eq!(
    //     poly,
    //     DenseUnivariatePolynomial::new(vec![Fq::from(0), Fq::from(0), Fq::from(1)])
    // );
}
