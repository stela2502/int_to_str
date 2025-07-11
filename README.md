[![Rust](https://github.com/stela2502/int_to_str/actions/workflows/rust.yml/badge.svg?branch=main)](https://github.com/stela2502/int_to_str/actions/workflows/rust.yml)

# int_to_str is a Rust DNA to 2bit encoding tool

It is able to en- and de-code DNA as 2bit integer. It's implementation is likely crappy without end, but it seams to be working fine.

# Usage


Add this to yout Cargo.toml
```
[dependencies]
int_to_str = { git="https://github.com/stela2502/int_to_str.git"}
```

In your code you can encode a DNA string as 2bit using

```

use int_to_str::int_to_str::IntToStr;

let tool = IntToStr::new(b"ATGACTCTCAGCATGGAAGGACAGCAGAGACCAAGAGATCCTCCCACAGGGACACTACCTCTGGGCCTGGGATAC");
```

You can convert data to any uint integer starting from the first bp or you can get a slice of the sequence (as IntToSeq) using the slice( start:usize, end:usize) function.

# Binary

The library also codes for a small binary:

```
int_to_str 
Usage: target/release/int_to_str <sequence or integer>
```

It either converts an integer (max u128) to DNS seq like this

```
int_to_str 12343
Integer input: 12343
→ Sequence: TCTAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
```

or converts a sequence to an integer (u128):
```

int_to_str TCTAAAT
Sequence input: TCTAAAT
→ Sequence: 12343

```

You see you need to manually make sure you only use the significant bits of your int as DNA seq - the binary does not help you there ;-)
