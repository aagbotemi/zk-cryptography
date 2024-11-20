use ark_ff::PrimeField;
use std::collections::HashMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Column {
    LEFT,
    RIGHT,
    OUTPUT,
}

impl From<usize> for Column {
    fn from(value: usize) -> Self {
        match value {
            1 => Column::LEFT,
            2 => Column::RIGHT,
            3 => Column::OUTPUT,
            _ => panic!("wrong column"),
        }
    }
}

#[derive(Debug, Clone)]
pub struct Cell {
    pub column: Column,
    pub row: usize,
}

impl Cell {
    pub fn label<F: PrimeField>(&self, group_order: u64) -> F {
        roots_of_unity::<F>(group_order)[self.row]
            * match self.column {
                Column::LEFT => F::one(),
                Column::RIGHT => F::from(2u8),
                Column::OUTPUT => F::from(3u8),
            }
    }
}
pub fn root_of_unity<F: PrimeField>(group_order: u64) -> F {
    F::get_root_of_unity(group_order).unwrap()
}

pub fn roots_of_unity<F: PrimeField>(group_order: u64) -> Vec<F> {
    let mut res = vec![F::from(1u32)];
    let generator: F = root_of_unity(group_order);
    for _ in 1..group_order {
        res.push(res[res.len() - 1] * generator);
    }
    res
}

pub fn merge_maps<F: PrimeField>(
    map1: &HashMap<Option<String>, F>,
    map2: &HashMap<Option<String>, F>,
) -> HashMap<Option<String>, F> {
    let mut merged = HashMap::new();
    for (key, val) in map1.iter().chain(map2.iter()) {
        *merged.entry(key.clone()).or_insert(F::zero()) += val;
    }
    merged
}

pub fn multiply_maps<F: PrimeField>(
    map1: &HashMap<Option<String>, F>,
    map2: &HashMap<Option<String>, F>,
) -> HashMap<Option<String>, F> {
    let mut result = HashMap::new();
    for (k1, v1) in map1.iter() {
        for (k2, v2) in map2.iter() {
            let product_key = get_product_key(k1.clone(), k2.clone());
            *result.entry(product_key).or_insert(F::zero()) += *v1 * v2;
        }
    }
    result
}

pub fn get_product_key(key1: Option<String>, key2: Option<String>) -> Option<String> {
    match (key1, key2) {
        (Some(k1), Some(k2)) => {
            let members = {
                let mut members = Vec::new();
                members.extend(k1.split('*'));
                members.extend(k2.split('*'));
                members.sort();
                members
            };
            Some(
                members
                    .into_iter()
                    .filter(|x| !x.is_empty())
                    .collect::<Vec<&str>>()
                    .join("*"),
            )
        }
        (Some(k1), None) => Some(k1),
        (None, Some(k2)) => Some(k2),
        (None, None) => None,
    }
}

pub fn is_valid_variable_name(name: &str) -> bool {
    !name.is_empty()
        && name.chars().all(char::is_alphanumeric)
        && !name.chars().next().unwrap().is_numeric()
}

pub fn evaluate<F: PrimeField>(exprs: &Vec<&str>) -> HashMap<Option<String>, F> {
    evaluate_inner(exprs, false)
}

pub fn evaluate_inner<F: PrimeField>(
    exprs: &Vec<&str>,
    first_is_negative: bool,
) -> HashMap<Option<String>, F> {
    match exprs.iter().any(|&x| x == "+") {
        true => {
            let idx = exprs.iter().position(|&x| x == "+").unwrap();
            let l = evaluate_inner(&exprs[..idx].to_vec(), first_is_negative);
            let r = evaluate_inner(&exprs[idx + 1..].to_vec(), false);
            return merge_maps(&l, &r);
        }
        false => match exprs.iter().any(|&x| x == "-") {
            true => {
                let idx = exprs.iter().position(|&x| x == "-").unwrap();
                let l = evaluate_inner(&exprs[..idx].to_vec(), first_is_negative);
                let r = evaluate_inner(&exprs[idx + 1..].to_vec(), true);
                return merge_maps(&l, &r);
            }
            false => match exprs.iter().any(|&x| x == "*") {
                true => {
                    let idx = exprs.iter().position(|&x| x == "*").unwrap();
                    let l = evaluate_inner(&exprs[..idx].to_vec(), first_is_negative);
                    let r = evaluate_inner(&exprs[idx + 1..].to_vec(), first_is_negative);
                    return multiply_maps(&l, &r);
                }
                false => {
                    if exprs.len() > 1 {
                        panic!("No ops, expected sub-expr to be a unit: {:?}", exprs[1]);
                    } else if exprs[0].starts_with('-') {
                        return evaluate_inner(&vec![&exprs[0][1..]], !first_is_negative);
                    } else if exprs[0].parse::<i128>().is_ok() {
                        let value = {
                            if first_is_negative {
                                F::from(exprs[0].parse::<i128>().unwrap().abs() as u128).neg()
                            } else {
                                F::from(exprs[0].parse::<i128>().unwrap() as u128)
                            }
                        };
                        let mut result = HashMap::new();
                        result.insert(None, value);
                        return result;
                    } else if is_valid_variable_name(exprs[0]) {
                        let mut result = HashMap::new();
                        let value = if first_is_negative {
                            F::one().neg()
                        } else {
                            F::one()
                        };
                        result.insert(Some(exprs[0].to_string()), value);
                        return result;
                    } else {
                        println!("exprs:{:?}", exprs);
                        panic!("ok wtf is {}", exprs[0]);
                    }
                }
            },
        },
    }
}
