use std::collections::BTreeMap;
use std::fmt;
//use crate::errors::SeqError;
//use crate::traits::BinaryMatcher;

// logics copied from https://github.com/COMBINE-lab/kmers/
pub type Base = u8;
pub const A: Base = 0;
pub const C: Base = 1;
pub const G: Base = 2;
pub const T: Base = 3;

#[derive(Default,Debug,PartialEq)]
pub struct IntToStr {
	long_term_storage: Vec::<u8>, //this will never be shifted nor poped
	storage: Vec::<u8>,  // the initial utf8 encoded data
	pub u8_encoded: Vec::<u8>, // the 2bit encoded array (4 times compressed)
	pub lost:usize, // how many times have I lost 4bp?
	pub shifted:usize, // how many times have we shifted our initial sequence?
	pub kmer_size:usize,
	step_size:usize,
	checker:BTreeMap::<u8, usize>,
	mask:u64, //a mask to fill matching sequences to match the index's kmer_len
	current_position:usize,

}

use std::mem::size_of;


macro_rules! impl_to_le_bytes {
    ($($t:ty),*) => {
        $(impl ToLeBytes for $t {
            fn to_le_byte_vec(&self) -> Vec<u8> {
                self.to_le_bytes().to_vec()
            }
        })*
    };
}

impl_to_le_bytes!(u8, u16, u32, u64, u128);

impl fmt::Display for IntToStr {

    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {

    	// Helper: format u64 or usize into grouped binary string
        fn to_bin_grouped<T: Into<u64>>(num: T, bits: usize) -> String {
            let raw = format!("{:0width$b}", num.into(), width = bits);
            raw.chars()
                .collect::<Vec<_>>()
                .chunks(4)
                .map(|chunk| chunk.iter().collect::<String>())
                .collect::<Vec<_>>()
                .join("_")
        }

        // Helper: format u8 vector as binary
        fn vec_to_bin(vec: &Vec<u8>) -> Vec<String> {
            vec.iter()
                .map(|b| format!("0b{}", to_bin_grouped(*b, 8)))
                .collect()
        }

        let dna_str = match std::str::from_utf8(&self.long_term_storage) {
            Ok(s) => s,
            Err(_) => "<invalid UTF-8>",
        };

        writeln!(f, "IntToStr {{")?;
        writeln!(f, "  long_term_storage: {:?},", self.long_term_storage)?;
 		writeln!(f, "  long_term_storage (DNA): \"{}\",", dna_str)?;
        writeln!(f, "  storage: {:?},", self.storage)?;
        writeln!(f, "  u8_encoded: {},", vec_to_bin(&self.u8_encoded).join(", "))?;
        writeln!(f, "  lost: {},", self.lost)?;
        writeln!(f, "  shifted: {},", self.shifted)?;
        writeln!(f, "  kmer_size: {},", self.kmer_size)?;
        writeln!(f, "  step_size: {},", self.step_size)?;
        writeln!(f, "  checker: {:?},", self.checker)?;
        writeln!(f, "  mask: 0x{},", to_bin_grouped(self.mask, 64))?; // hex for clarity
        writeln!(f, "  current_position: {}", self.current_position)?;
        write!(f, "}}")
    }
}

// Implement the Index trait for MyClass
use std::ops::Index;

impl Index<usize> for IntToStr {
    type Output = u8;

    fn index(&self, index: usize) -> &Self::Output {
        self.u8_encoded.get(index).unwrap_or(&0_u8)
    }
}


impl PartialEq<Vec<u8>> for IntToStr {
    fn eq(&self, other: &Vec<u8>) -> bool {
        &self.u8_encoded == other
    }
}

/*
impl Iterator for IntToStr {
    type Item = (u16, u64);

    fn next(&mut self) -> Option<Self::Item> {

        let result = self.seq_at_position(self.current_position);

        if let Some((cell_id, second_seq)) = result {
            // Advance to the next position
            self.current_position += 8;
            Some((cell_id, second_seq))
        } else {
            // If the shift is 7, stop the iteration - we reached the end
            if (self.current_position % 8) + self.step_size >= 8  {
                return None
            }
            // Otherwise, move to the next position
            self.current_position = (self.current_position % 8) + self.step_size;

            // Get the result at the updated position
            let result = self.seq_at_position(self.current_position);

	        if let Some((cell_id, second_seq)) = result {
	            // Advance to the next position
	            self.current_position += 8;
	            Some((cell_id, second_seq))
	        } else {
	        	// reset the iterator for the next use
	        	self.current_position = 0;
	            None
	        }
        }
    }
}
*/
/// Here I have my accessorie functions that more or less any of the classes would have use of.
/// I possibly learn a better way to have them...
impl IntToStr {

	pub fn new(seq: &[u8] ) -> Self {
		// 4 of the array u8 fit into one result u8
		//eprintln!("Somtimes I die?! -> processed seq: {:?}", seq);
		fn encode_binary( c: u8) -> Result<Base, String> {
	    	// might have to play some tricks for lookup in a const
	    	// array at some point
	    	match c {
	        	b'A' | b'a' => Ok(A),
	        	b'C' | b'c' => Ok(C),
		        b'G' | b'g' => Ok(G),
	        	b'T' | b't' => Ok(T),
	        	b'N' | b'n' => Ok(A), // this is necessary as we can not even load a N containing sequence
		        _ => Err("cannot encode {c} into 2 bit encoding".to_string()),
	    	}
		}

		let storage:Vec::<u8> = seq.to_vec();
		let long_term_storage = seq.to_vec();

		let num_bytes =(storage.len()+3)/4;
		let mut u8_encoded = Vec::<u8>::with_capacity( num_bytes );
		for id in 0..num_bytes {
			u8_encoded.push( 0_u8 );
			for add in (0..4).rev(){
				let current_utf8_id = id*4 + add;
				u8_encoded[id] <<=2;
				if seq.len() > current_utf8_id {
					u8_encoded[id] |= match encode_binary(seq[current_utf8_id]){
						Ok(bits) => bits,
						Err(e) => panic!("Char at seq[{current_utf8_id}] ({}) is not a supported nucleotoide",seq[current_utf8_id] )
					};
				} // this automatically "stuffs" the last byte with zeros's
			}
		}
		let lost = 0;
		let shifted = 0;
		let checker = BTreeMap::<u8, usize>::new();
		//let mask: u64 = (1 << (2 * kmer_size)) - 1;

		Self{
            long_term_storage,
			storage,
			u8_encoded,
			lost,
			shifted,
			kmer_size:16, //deprecated
			checker,
			mask: 0, // useless
			current_position:0,
			step_size:1,
		}
	}

	pub fn step_size( &mut self, size:usize ) {
		if size > 7 {
			eprintln!("step_size can not be larger than 7");
			self.step_size = 7;
		}else if size == 0 {
			eprintln!("step_size can not be zero");
			self.step_size = 1;
		}else {
			self.step_size = size;
		}
	}


    // This function retures a Option<( u16, SecondSeq )> - an enhanced u16 and an enhanced u64
    // or None if the u64 would only consist of less than 4 bytes as the remaining sequence is to short.
    pub fn seq_at_position(&self, start:usize ) -> Option<( u16, u64 )> {
    	// Ensure the start position is within bounds
	    if start >= self.storage.len() {
	        return None;
	    }
    	
    	// Determine the range of bytes to consider
    	let start_id = start / 4;
    	let stop:u8;
    	let shift_amount = (start % 4) * 2;

    	// keep the later check from getting negative
    	if self.storage.len() <= 8 + start {
    		return None
    	}

    	let stop_id = if start + 8 + self.kmer_size + shift_amount  < self.storage.len() {
    		stop = self.kmer_size as u8;
		    start_id + 2 + (self.kmer_size + 3 + shift_amount) / 4
		} else {
			//println!("{start} + 8 + {} +3 + {shift_amount}  < {}",self.kmer_size, self.storage.len() );
			stop = (self.storage.len() - 8 - start) as u8;
			//println!( "{stop} = ({} - 8 - {start}) as u8;", self.storage.len());
		    self.u8_encoded.len()
		};

		let min_length = self.kmer_size.min(20);
		
    	if (stop as usize) < min_length  {
    		//println!("seq_at_position - to small u64: {stop} < {min_length}");
    		return None
    	}


    	//println!("I got start_id {start_id} and stop_id {stop_id}");
    	// Create a u128 from the selected byte range
	    let mut u128_value = 0u128;
	    for byte in self.u8_encoded[start_id..stop_id].iter().rev() {
	        u128_value <<= 8;
	        u128_value |= *byte as u128;
	    }

	    //println!( "and I got the the u128 {u128_value:b}");
	    // Shift the u128 value to the right position
	    
	    let shifted_value = u128_value >> shift_amount;
	    
	    //println!("I got the shift amount {shift_amount} [bits] and the u64 len {stop}");

	    // Extract u16 and u64 values
	    let u16_value = (shifted_value & u16::MAX as u128) as u16;
	    let u64_value = (shifted_value >> 16) as u64;

	    //println!("From that u128 I created \nu16: {u16_value:b} and\nu64 {u64_value:b}");

	    // Return the result
    	Some( (u16_value, u64_value ) )
    }

    pub fn is_empty(&self) -> bool{
    	self.u8_encoded.is_empty()
    }

	pub fn len(&self) -> usize{
		self.u8_encoded.len()
	}
	
	fn reverse_bits_in_byte( b: u8) -> u8 {
		let mut b = b;
	    b = (b >> 4) | (b << 4);
	    b = ((b & 0b1100_1100) >> 2) | ((b & 0b0011_0011) << 2);
	    b = ((b & 0b1010_1010) >> 1) | ((b & 0b0101_0101) << 1);
	    b
	}

	///helper function for crating the unsigned intergers from the internal u8 array.
    fn into_uint<T, const N: usize>(
	    encoded: &[u8],
	    reverse_bits: bool,
	) -> T
	where
	    T: Default
	        + From<u8>
	        + std::ops::Shl<u8, Output = T>
	        + std::ops::BitOr<Output = T>,
	{
	    let mut ret = T::default();

	    for i in (0..N).rev() {
	        let mut byte = encoded.get(i).copied().unwrap_or(0);
	        if reverse_bits {
	            byte = Self::reverse_bits_in_byte(byte);
	        }

	        ret = (ret << 8) | T::from(byte);
	    }

	    ret
	}
	/// this is needed for the initial mappers
	/// takes the self encoded seqences and converts the last entries into one u16
	pub fn into_u16(&self ) -> u16{
		Self::into_uint::<u16, 2>(&self.u8_encoded, false)
	}	
	/// this is needed for completing the set
	pub fn into_u32(&self ) -> u32{
		Self::into_uint::<u32, 4>(&self.u8_encoded, false)
	}
	/// needed for the secondary mappings
    /// takes the UTF8 encoded sequence and encodes the first 32 into a u64 
	pub fn into_u64(&self ) -> u64{
		Self::into_uint::<u64, 8>(&self.u8_encoded, false)
	}
	/// this is needed for completing the set
	pub fn into_u128(&self ) -> u128{
		Self::into_uint::<u128, 16>(&self.u8_encoded, false)
	}


	/// reverse into_uint - create the internal structure from the numbers
	fn from_uint<T, const N: usize>(mut value: T, reverse_bits: bool) -> Self
	where
	    T: Copy + ToLeBytes
	        + std::ops::Shr<u8, Output = T>
	        + std::ops::BitAnd<Output = T>
	        + From<u8>
	        + Into<u128>, // For size handling
	{
		let bytes = value.to_le_byte_vec();

		// Ensure we have exactly N bytes, padding with 0s if needed
	    let mut buf = vec![0u8; N];
	    for (i, b) in bytes.iter().enumerate().take(N) {
	        buf[i] = *b;
	    }

	    if reverse_bits {
	        for b in &mut buf {
	            *b = Self::reverse_bits_in_byte(*b);
	        }
	    }
		let utf8 = Self::u8_array_to_str(&bytes);

	    IntToStr{
	    	u8_encoded: buf,
	    	storage: utf8.clone().into_bytes(), 
			long_term_storage: utf8.into_bytes(),
			lost: 0,
			shifted: 0,
			kmer_size: 16,
			step_size: 3,
			checker: BTreeMap::new(),
			mask: 0_u64,
			current_position: 0,
	    }
	}
	
	pub fn from_u16(val: u16) -> Self {
	    Self::from_uint::<u16, 2>(val, false)
	}

	pub fn from_u32(val: u32) -> Self {
	    Self::from_uint::<u32, 4>(val, false)
	}

	pub fn from_u64(val: u64) -> Self {
	    Self::from_uint::<u64, 8>(val, false)
	}

	pub fn from_u128(val: u128) -> Self {
	    Self::from_uint::<u128, 16>(val, false)
	}

    /// needed for the secundary mappings
    /// takes the UTF8 encoded sequence and encodes the first kmer_len into a u64 
    pub fn into_u64_nbp(&self, kmer_size:usize ) -> u64{

        let mut ret: u64 = 0;
        let target = kmer_size / 4;
        let residual = kmer_size % 4;

        if residual != 0 {
        	println!("residual == {residual} - fill up with zeros!");
        	if self.u8_encoded.len() <= target {
        		ret <<= 8;
        	}
        	else {
        		//println!( "have {} and want {}", self.u8_encoded.len() , target);
        		let mut loc = self.u8_encoded[target];
	        	//println!("\nGet the loc from self.u8_encoded[{}].clone()", target );
	        	//println!("I got this u64 for my u8: {loc:b}, {:b}", self.u8_encoded[target]);
	        	match residual{
	        		0 => (),
	        		3 => loc &= 0b00000011,
	        		2 => loc &= 0b00001111,
	        		1 => loc &= 0b00111111,
	        		_ => panic!("what are you trying to atchieve here? {residual}"),
	        	};
	        	ret = (ret << 8) | loc as u64;
        	}
        	
        }

        for i in (0..target).rev() {
            ret = (ret << 8) | self.u8_encoded.get( i ).copied().unwrap_or(0) as u64;
        }

        

        ret
    }




	/// drop the last n u8 2 bit encoded sequences (n=1 => 4bp dropped from encoded)
	/// this can be reset using the self.reset() function.
	pub fn dropped_n (&mut self, n:usize ) -> Option<()>{
		//let mut  removed:String = "".to_string();
		//self.print();
		if self.u8_encoded.len() < 2{
			return None;
		}
		for _i in 0..n{
			let _a =self.u8_encoded.remove(0);
			//removed.clear();
			//self.u8_to_str( 4, &a, &mut removed);
			//eprintln!("dropped_n {i} has dropped {:?}",removed);
		}
		self.lost +=n;
		//self.print();
		Some ( () )
	}

	/// shift the underlying sequence by 1 (frameshift)
	/// can be regenrated using deep_refresh().
	pub fn shift(&mut self )-> Option<()> {
		panic!("This functionallyit has been removed!");
	}


	/// Convert any integer value obtained from any of the self.into_u16.to_le_bytes() like functions
	/// and adds them at the end of the mutable data string.
	pub fn to_string (&self, bases:usize )->String{
		let mut i = 0;
		let mut data = String::default();
	    
		for u8_4bp in self.u8_encoded.iter(){
			i += 4;
			Self::u8_to_str(  u8_4bp, &mut data );
			if i >= bases {				
				break;
			}
		}
		data.truncate(bases);
		data
	}

	pub fn u8_to_str ( u8_rep:&u8,  data:&mut String ){

		let mut loc:u8 = *u8_rep;
		//println!("converting u8 {loc:b} to string with {kmer_size} bp.");
		for _i in (0..4).rev(){
			let ch = match loc & 0b11{
	            0 => "A", //0b00
	            1 => "C", // 0b01
	            2 => "G", // 0b10
	            3 => "T", // 0b11
	            _ => "N",
	       	};
	       	*data += ch;
	       	loc >>= 2;
	       	//println!("{ch} and loc {loc:b}");
		}

		//println!("\nMakes sense? {:?}", data);
	}

	pub fn u8_array_to_str ( u8_rep:&[u8] ) ->String{

		let mut i = 0;
		let mut data = "".to_string();

		for u8_4bp in u8_rep {
			Self::u8_to_str( u8_4bp, &mut data );
		}
		data
	}

	/// This will mask the providied u64 with the internal mask definitivels 'killing' the last bp.
	/// This will convert the masked bp into 'A' or better to say #b00 entries.
	pub fn mask_u64( &self, seq:&u64 ) ->u64{
		//let filled_sed = seq | (!self.mask & (0b11 << (2 * self.kmer_size )));
		//println!("I have the mask {:b}", self.mask);
        //seq & self.mask
        *seq
	}

	/// This will mask the last bp by 'A' or #b00 and only return bp length usabel info.
	pub fn mask_u64_to ( &self, seq:&u64, bp:usize) -> u64  {
		*seq
		/*if bp >= 32{
			return *seq;
		}
		let mask = !0u64 >> ( (32 - bp )  * 2) as u32;
		seq & mask*/
	}


	pub fn print(&self) {
		println!(">seq (n={} [{}*4 + {}])\n{}", self.storage.len(),  self.storage.len()/4,self.storage.len() %4, std::str::from_utf8( &self.storage.clone() ).expect("Invalid UTF-8"));
		println!("I have 'lost' {} u8's from my internal encoded", self.lost );
		print!("[ ");
		for en in &self.u8_encoded{
			print!(", {en:b}");
		}
		println!("]");
		println!("human readable: {:?}", self.u8_encoded );
		println!("binary vecor length {}", self.len());
	}

}


/*---- int_to_str::tests::test_uXx_roundtrip_aaaaaaacaaagaaat stdout ----
==> Testing input: AAAAAAACAAAGAAAT
Original encoded: [0, 1, 2, 3]
Encoded u16:                   1000000000000000
Encoded u32:   11000000 01000000 10000000 00000000
Encoded u64:   11000000 01000000 10000000 00000000
Encoded u128:  11000000 01000000 10000000 00000000
Decoded u16:   00000000 00000001
Decoded u32:   00000000 00000001 00000010 00000011
Decoded u64:   00000000 00000001 00000010 00000011, 00000000, 00000000, 00000000, 00000000
Decoded u128:  00000000 00000001 00000010 00000011, 00000000, 00000000, 00000000, 00000000, 00000000, 00000000, 00000000, 00000000, 00000000, 00000000, 00000000, 00000000
*/


#[cfg(test)]
mod tests {
    use super::*;
    use crate::IntToStr;



	fn format_bytes_binary(bytes: &[u8]) -> String {
	    bytes.iter()
	        .map(|b| format!("{:08b}", b))
	        .collect::<Vec<_>>()
	        .join(", ")
	}
	fn roundtrip_test(input: &str) {
	    
	    let used_str = if input.len() < 64 {
	    	let padding = "A".repeat(64 - input.len());
        	format!("{input}{padding}")
	    }else {
	    	input.to_string()
	    };
		println!("==> Testing input: {used_str}");
	    let encoder = IntToStr::new(used_str.as_bytes());
	    let original = encoder.u8_encoded.clone();
	    println!("Original encoded: {:?}", original);

	    // Encode
	    let val16 = encoder.into_u16();
	    let val32 = encoder.into_u32();
	    let val64 = encoder.into_u64();
	    let val128 = encoder.into_u128();

	    println!("Encoded u16:   {:b}", val16);
	    println!("Encoded u32:   {:b}", val32);
	    println!("Encoded u64:   {:b}", val64);
	    println!("Encoded u128:  {:b}", val128);

	    //panic!("Will die anyhow");
	    // Decode
	    let decoded16 = IntToStr::from_u16(val16);
	    let decoded32 = IntToStr::from_u32(val32);
	    let decoded64 = IntToStr::from_u64(val64);
	    let decoded128 = IntToStr::from_u128(val128);

	    println!("Decoded u16:   {}", format_bytes_binary(&decoded16.u8_encoded));
	    println!("Decoded u32:   {}", format_bytes_binary(&decoded32.u8_encoded));
	    println!("Decoded u64:   {}", format_bytes_binary(&decoded64.u8_encoded));
	    println!("Decoded u128:  {}", format_bytes_binary(&decoded128.u8_encoded));

	    // Assertions
	    assert_eq!(&original[..2], &decoded16.u8_encoded[..2], "u16 mismatch");
	    assert_eq!(&original[..4], &decoded32.u8_encoded[..4], "u32 mismatch");
	    assert_eq!(&original[..8], &decoded64.u8_encoded[..8], "u64 mismatch");
	    assert_eq!(&original[..16], &decoded128.u8_encoded[..16], "u128 mismatch");

	    //assert_eq!(&encoder.storage[..8], &decoded16.storage[..8], "u16 mismatch");
	    //assert_eq!(&encoder.storage[..16], &decoded32.storage[..16], "u32 mismatch");
	    //assert_eq!(&encoder.storage[..32], &decoded64.storage[..32], "u64 mismatch");
	    assert_eq!( &decoded128.to_string( input.len() ), input );

	}
	#[test]
	fn test_encoding(){
		let encoder = IntToStr::new(b"ACGTTC");
		assert_eq!( encoder.into_u16(), 0b011111100100, "encoded {:b} is not the expected {:b}",
		 encoder.into_u16(), 0b011111100100 );
		assert_eq!( &encoder.to_string( 5 ), "ACGTT" );
		//this will only overreach by max 3 bp of A's (0)
		assert_eq!( &encoder.to_string( 15 ), "ACGTTCAA" );
	}

	#[test]
	fn test_from_uint(){
		// the from_uint function is the 'mother' of all fromn16 to from_u128.
		// So testing one should also test all others.
		let encoder = IntToStr::from_u16( 0b1111100100 as u16 );
		let exp = IntToStr::new(b"ACGTT");
		assert_eq!( encoder.u8_encoded, exp.u8_encoded );
	}

    #[test]
    fn test_uXx_roundtrip_aaaaaaacaaagaaat() {
        roundtrip_test("AAAAAAACAAAGAAAT");
    }

    #[test]
    fn test_uXx_roundtrip_aaaaaacaagaataaa() {
        roundtrip_test("AAAAAACAAGAATAAA");
    }
}