pub mod protocol;
pub mod succint_protocol;
pub mod utils;

pub fn play_game(n: u32, print: bool) {
    let result = fizz_buzz(n);
    if print {
        // println!("{result}");
    }
}

pub fn fizz_buzz(n: u32) -> String {
    match (n % 3, n % 5) {
        (0, 0) => "FizzBuzz".to_string(),
        (0, _) => "Fizz".to_string(),
        (_, 0) => "Buzz".to_string(),
        (_, _) => n.to_string(),
    }
}
