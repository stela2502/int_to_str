// IntToString.rs

// Define your library code here

// Include a module for tests
#[cfg(FALSE)]

#[cfg(test)]
mod tests {

   use int_to_str::int_to_str::IntToStr;
   //use int_to_str::int_to_str::IntToStr;


   #[test]
    fn test_seq_at_position() {
      let tool = IntToStr::new(b"ATGACTCTCAGCATGGAAGGACAGCAGAGACCAAGAGATCCTCCCACAGGGACACTACCTCTGGGCCTGGGATAC");
      if let Ok( ( cellid, second_seq) ) = tool.seq_at_position(0){
         panic!( "{:b} - {:b}",cellid, second_seq);
         let seq_u16 = tool.to_string(8);
         //println!("The sequenc I got: {cellid:b} should be ATGACTCT");
         assert_eq!( seq_u16, "ATGACTCT".to_string() );
         let seq_u64 = tool.to_string( 32, &second_seq );
         assert_eq!( seq_u64, "CAGCATGGAAGGACAGCAGAGACCAAGAGATC".to_string() );
      }else {
         panic!("seq_at_position did not return a value")
      }
   }

   #[test]
    fn test_seq_at_position_smaller_u64() {
      let tool = IntToStr::new(b"ATGACTCTCAGCATGGAAGGACAGCAGAGACCAAGAGA");
      if let Ok( ( cellid, second_seq) ) = tool.seq_at_position(0){
         let seq_u16 = tool.u64_to_string( 8, &(cellid as u64));
         //println!("The sequenc I got: {cellid:b} should be ATGACTCT");
         assert_eq!( seq_u16, "ATGACTCT".to_string() );
         let seq_u64 = tool.u64_to_string(30 , &second_seq );
         assert_eq!( seq_u64, "CAGCATGGAAGGACAGCAGAGACCAAGAGA".to_string() );
      }else {
         panic!("seq_at_position did not return a value")
      }
   }

   #[test]
    fn test_seq_at_position_start_3() {
      let tool = IntToStr::new(b"ATGACTCTCAGCATGGAAGGACAGCAGAGACCAAGAGATCCTCCCACAGGGACACTACCTCTGGGCCTGGGATAC");
      if let Ok( ( cellid, second_seq) ) = tool.seq_at_position(3){
         let seq_u16 = tool.u64_to_string( 8, &(cellid as u64));
         //println!("The sequenc I got: {cellid:b} should be ATGACTCT");
         assert_eq!( seq_u16, "ACTCTCAG".to_string() );
         let seq_u64 = tool.u64_to_string( 32, &second_seq );
         assert_eq!( seq_u64, "CATGGAAGGACAGCAGAGACCAAGAGATCCTC".to_string() );
      }else {
         panic!("seq_at_position did not return a value")
      }
   }

   #[test]
    fn test_seq_at_position_out_of_range() {
      let tool = IntToStr::new(b"ATGACTCTCAGCATGGAAGGACAGCAGAGACCAAGAGATCCTCCCACAGGGACACTACCTCTGGGCCTGGGATAC");
      match tool.seq_at_position(70){
         Some( ( cellid, second_seq) ) => panic!("expected None for an out of range id! And got {cellid} and {second_seq}"),
         None => assert!(true),
      };
   }

    #[test]
    fn test_u64_to_str(){

        let tool = IntToStr::new( b"CGATATT");

        let num:u64 = tool.into_u64();
        println!("I have this number for the sting 'CGATATT' {num}");

        //let num:u64 = 15561;

        //println!("This is the number I want to convert to a Sequence: {num:b}");
        // 11110011001001
        // T T A T A C G 
        //let tool = IntToStr::new();
        let mut data:String = "".to_string();
        tool.to_string( 7, &mut data );
        assert_eq!( data, "CGATATT".to_string() )
    } 

     #[test]
    fn check_conversion_4bp() {

     let seq = b"AGGC";
     //         C G G A
     //         01101000
     let tool = IntToStr::new( seq );

     assert_eq!( tool.len(),  1 ); 
     //panic!("{:b}", binary[0] );
     //                          G G C 
     assert_eq!( tool, vec![ 0b1101000 ]);

     let mut decoded:String = "".to_string();

     tool.to_string( 15, &mut decoded );
     assert_eq!( decoded, "AGGC" );      
                            
    }


     #[test]
    fn check_conversion_15bp() {
        //          0000111122223333   
     let seq = b"AGGCTTGATAGCGAG";
     let tool = IntToStr::new(seq);

     assert_eq!( tool.len(),  4 );

     //panic!("{:b} {:b} {:b} {:b}", binary[0], binary[1], binary[2], binary[3] );
     //                                                A G C A
     //                                              0b00100100  
     //                          G G C     T T G A     T A G C     G A G

     //panic!("the binaries I get: {:b} {:b} {:b} {:b} ", binary[0], binary[1], binary[2], binary[3]);
     //assert_eq!( binary, vec![ 0b1101000, 0b101111, 0b1100011, 0b100010 ]);
     //tool.print();
     assert_eq!( tool, vec![ 0b1101000, 0b101111, 0b1100011, 0b100010  ]);
     let mut decoded:String = "".to_string();

     tool.to_string(15, &mut decoded );
     assert_eq!( decoded, "AGGCTTGATAGCGAG" );
    }

     #[test]
    fn check_conversion_1bp() {

     let seq = b"C";
     let tool = IntToStr::new( seq );

     assert_eq!( tool.len(),  1 ); 
     //                                                A G C A
     //                                              0b00100100  
     //                          G G C     T T G A     T A G C     G A G
     assert_eq!( tool, vec![ 0b1 ]);

     let mut decoded:String = "".to_string();

     tool.to_string(1, &mut decoded );
     assert_eq!( decoded, "C" );
    }

    #[test]
    fn check_conversion_one_a() {

     let seq = b"A";
     let tool = IntToStr::new( seq );

     assert_eq!( tool.len(),  1 ); 
     //                                                A G C A
     //                                              0b00100100  
     //                          G G C     T T G A     T A G C     G A G
     assert_eq!( tool, vec![ 0b0 ]);

     let mut decoded:String = "".to_string();

     tool.to_string(1, &mut decoded );
     assert_eq!( decoded, "A" );
    }


    #[test]
    fn check_conversion_4_from_15bp() {
     //          ----    ----
     let seq = b"AGGCCTGTATGA";
     let tool = IntToStr::new( seq );

     assert_eq!( tool.len(),  3 ); 

     //                     A G T A 
     tool.print();
     assert_eq!( tool[2], 0b101100 );

     let mut decoded:String = "".to_string();

        tool.to_string(4, &mut decoded );
        assert_eq!( decoded, "AGGC" );
        decoded.clear();
        tool.to_string(3, &mut decoded );
     assert_eq!( decoded, "AGG" );                
    }

    #[test]
    fn check_longer_string() {

     let seq = b"CTGGAAGCGCTGGGCTCCCGGCTGCATTGGGCTGGTCCGTGGGTC";
     let tool = IntToStr::new(seq );

     assert_eq!( tool.len(), 12 ); 

     //                       CT G G
     //                       G G T C
     tool.print();
     assert_eq!( tool[11], 0b1  );

     let mut decoded:String = "".to_string();

        tool.to_string(4,  &mut decoded );
        assert_eq!( decoded, "CTGG" );      
        decoded.clear();
        tool.to_string(tool.len()*4 -3, &mut decoded );
        assert_eq!( decoded, "CTGGAAGCGCTGGGCTCCCGGCTGCATTGGGCTGGTCCGTGGGTC" );       

        decoded.clear();
        tool.to_string(tool.len()*4 -6,  &mut decoded );
        assert_eq!( decoded, "CTGGAAGCGCTGGGCTCCCGGCTGCATTGGGCTGGTCCGTGG" );               
    }

    #[test]
    fn check_u64_decode() {

     let seq = b"CTGGAAAAGCTGGGCTCCCGGCTGCATTGGGCTGGTCCGTGGGTT";
     //let seq = b"CTGG";
     let mut tool = IntToStr::new(seq );
     tool.print();

     let two_bp = tool.into_u64_nbp( 2 );
     println!("the obtained sequence from into_u64_nbp for the sequence CT : {two_bp:b}");
     assert_eq!( tool.into_u64_nbp( 2 ) , 0b1101 );

     let binary = tool.into_u16(  ).to_le_bytes();
     println!("The u16 split ínto two u8 (again): {binary:?}", );
     assert_eq!( binary[0] as u8, tool.u8_encoded[0] );
     assert_eq!( binary[1] as u8, tool.u8_encoded[1] );
     

     //println!("Binary GGTC {:b}", IntToStr::enc_2bit(b"GGTC".to_vec() ) );
     assert_eq!( tool.u8_encoded, vec![173_u8, 0_u8, 182_u8, 218_u8, 149_u8, 182_u8, 241_u8, 106_u8, 235_u8, 229_u8, 234_u8, 3_u8]);


     let mut decoded:String = "".to_string();
     tool.to_string( 32, &mut decoded);
     assert_eq!( &decoded, "CTGGAAAAGCTGGGCTCCCGGCTGCATTGGGC" );
   
     decoded.clear();
     IntToStr::u8_to_str(  &tool.u8_encoded[0], &mut decoded );
     assert_eq!( &decoded, "CTGG" );

     decoded.clear();
     IntToStr::u8_array_to_str( &tool.u8_encoded);
     decoded.truncate( 45 );

     assert_eq!( &decoded, "CTGGAAAAGCTGGGCTCCCGGCTGCATTGGGCTGGTCCGTGGGTT" );

     decoded.clear();
     tool.to_string( 2, &mut decoded);
     assert_eq!( decoded, "CT".to_string() );

     let _ = tool.regenerate();
     // println!("If this is not followed by a three columns table we have a pporoblem!");
     // while let Some(t) = tool.next(){
     //     println!("{} {} {}", t.0, t.1, t.2);
     // }
     // tool.regenerate();
     let left = tool.next();
     println!("left: {:?}, right: {:?}", left, Some( (173_u16, 55990_u64) ) );
     let that = left.unwrap();
     assert_eq!( that.0, 173_u16 );
     assert_eq!( that.1, 55990_u64 );            
    }

    #[test]
    fn test_str_to_u64() {
      let seq = "AAAAAAAAAAAAAAAA"; // 16 bp A' => 0
      let mut tool = IntToStr::new(b"AGCT" );
      let numeric = tool.str_to_u64( seq );
      assert_eq!( numeric, 0_u64 );

      let seq2 = "AAAAAAACAAAGAAAT";
      println!("restart: {}",seq2);
      let numeric2 = tool.str_to_u64( seq2 ); 
      println!("binary result: {numeric2:b}");
      let mut res= "".to_string();
      tool.u64_to_str( 16, &numeric2, &mut res );
      println!("And back to seq: {res}",  );
      assert_eq!( seq2, res);
      assert_eq!( numeric2, 3229630464_u64 );

    }

    #[test]
    fn check_mask_u64() {
       let seq1_u64 = 14086662093597932231_u64;
       let tool = IntToStr::new(b"TAGTGTCCTGTGACTTCACCTCAAGTTGTAAT" );
       assert_eq!( seq1_u64, tool.into_u64(), "correct u64" );

       let mut seq = String::from("");
       tool.u64_to_str( 32, &seq1_u64, &mut seq);
       // last position will always be an A - sorry!
       assert_eq!( seq, String::from("TAGTGTCCTGTGACTTCACCTCAAGTTGTAAT"));
       eprintln!("This should be the masked seqence: \n{seq}\nTAGTGTCCTGTGACTTCACCTCAAGTTGTAAT"); 
       //println!("This should be the masked seqence: \n{seq}\nTAGTGTCCTGTGACTTCACCTCAAGTTGTAAT");
       //assert_eq!( masked, tool.into_u16() as u64, "correct masked 8bp u16" );
   }

   #[test]
    fn check_next() {
      let mut tool = IntToStr::new(
         b"ATGACTCTCAGCATGGAAGGACAGCAGAGACCAAGAGATCCTCCCACAGGGACACTACCTCTGGGCCTGGGATAC");
//          TGACTCTC       AAGGACAGCAGAGACCAAGAGATCCTCCCACAGGGACACT
      let mut first = "".to_string();
      let mut second = "".to_string();
      let mut i = 0;
      if let Some(entries) = tool.next(){
         i+=1;
         first.clear();
         second.clear();
         tool.u16_to_str( 8, &entries.0 , &mut first );
         assert_eq!( first, "ATGACTCT".to_string() );
         tool.u64_to_str( 32, &entries.1, &mut second );
      }else {
         assert_eq!( 1,2, "there should be a first entry in the tool!")
      }

      if let Some(entries) = tool.next(){
         i+=1;
         first.clear();
         second.clear();
         tool.u16_to_str( 8, &entries.0, &mut first );
         assert_eq!( first, "CAGCATGG".to_string() );
         tool.u64_to_str( 32, &entries.1, &mut second );
         assert_eq!( second, "AAGGACAGCAGAGACCAAGAGATCCTCCCACA".to_string() );
      }else {
         assert_eq!( 1,2, "there should be a second entry in the tool!")
      }

      if let Some(entries) = tool.next(){
         i+=1;
         first.clear();
         second.clear();
         tool.u16_to_str( 8, &entries.0, &mut first );
         assert_eq!( first, "AAGGACAG".to_string() );
         tool.u64_to_str( 32, &entries.1, &mut second );
         assert_eq!( second, "CAGAGACCAAGAGATCCTCCCACAGGGACACT".to_string() );
      }else {
         assert_eq!( 1,2, "there should be a second entry in the tool!")
      }
      if let Some(entries) = tool.next(){
         i+=1;
         first.clear();
         second.clear();
         tool.u16_to_str( 8, &entries.0, &mut first );
         assert_eq!( first, "CAGAGACC".to_string() );
         tool.u64_to_str( 32, &entries.1, &mut second );
         assert_eq!( second, "AAGAGATCCTCCCACAGGGACACTACCTCTGG".to_string() );
      }else {
         assert_eq!( 1,2, "there should be a third entry in the tool!")
      }

      if let Some(entries) = tool.next(){
         i+=1;
         first.clear();
         second.clear();
         tool.u16_to_str( 8, &entries.0, &mut first );
         assert_eq!( first, "AAGAGATC".to_string() );
         tool.u64_to_str( 32, &entries.1, &mut second );
         assert_eq!( second, "CTCCCACAGGGACACTACCTCTGGGCCTGGGA".to_string() );
      }else {
         assert_eq!( 1,2, "there should be a fourth entry in the tool!")
      }

      if let Some(entries) = tool.next(){
         i+=1;
         first.clear();
         second.clear();
         tool.u16_to_str( 8, &entries.0, &mut first );
         assert_eq!( first, "CTCCCACA".to_string() );
         tool.u64_to_str( 32, &entries.1, &mut second );
         assert_eq!( second, "GGGACACTACCTCTGGGCCTGGGATACAAAAA".to_string() );
      }else {
         assert_eq!( 1,2, "there should be a fifth entry in the tool!")
      }

      if let Some(entries) = tool.next(){
         i+=1;
         first.clear();
         second.clear();
         tool.u16_to_str( 8, &entries.0, &mut first );
         assert_eq!( first, "TGACTCTC".to_string() );
         tool.u64_to_str( 32, &entries.1, &mut second );
         assert_eq!( second, "AGCATGGAAGGACAGCAGAGACCAAGAGATCC".to_string() );
      }else {
         assert_eq!( 1,2, "there should be a sixth entry in the tool!")
      }

      /*if let Some(entries) = tool.next(){
         i+=1;
         first.clear();
         second.clear();
         tool.u16_to_str( 8, &entries.0, &mut first );
         assert_eq!( first, "AGCATGGA".to_string() );
         tool.u64_to_str( 32, &entries.1.0, &mut second );
         assert_eq!( second, "GCCTGGGATACAAAAAAAAAAAAAAAAAAAAA".to_string() );
         assert_eq!( entries.1.1, 11, "the expected max entry length" );
      }else {
         assert_eq!( 1,2, "there should be a sixth entry in the tool!")
      }*/

      while let Some(_entries) = tool.next(){
         i+=1;
      }
      assert_eq!( i,48, "A total of 54 fragments!")
   }

   #[test]
   fn test_antibody_tag(){
      let mut tool = IntToStr::new(
         b"CGAGAATTCCGATGCGCGTGTTAAGTATATAGGTTG");
      let mut first = "".to_string();
      let mut second = "".to_string();
      let mut i = 0;
      if let Some(entries) = tool.next(){
         i+=1;
         first.clear();
         second.clear();
         tool.u16_to_str( 8, &entries.0 , &mut first );
         assert_eq!( first, "CGAGAATT".to_string() );
         tool.u64_to_str( 28 , &entries.1, &mut second );
         assert_eq!( second, "CCGATGCGCGTGTTAAGTATATAGGTTG".to_string() );
      }else {
         assert_eq!( 1,2, "there should be a first entry in the tool!")
      }
      if let Some(entries) = tool.next(){
         i+=1;
         first.clear();
         second.clear();
         tool.u16_to_str( 8, &entries.0 , &mut first );
         assert_eq!( first, "CCGATGCG".to_string() );
         tool.u64_to_str( 20 , &entries.1, &mut second );
         assert_eq!( second, "CGTGTTAAGTATATAGGTTG".to_string() );
      }else {
         assert_eq!( 1,2, "there should be a second entry in the tool!")
      }
      if let Some(entries) = tool.next(){
         i+=1;
         first.clear();
         second.clear();
         tool.u16_to_str( 8, &entries.0 , &mut first );
         assert_eq!( first, "GAGAATTC".to_string() );
         tool.u64_to_str( 27, &entries.1, &mut second );
         assert_eq!( second, "CGATGCGCGTGTTAAGTATATAGGTTG".to_string() );
      }else {
         assert_eq!( 1,2, "there should be a third entry in the tool!")
      }
      while let Some(_entries) = tool.next(){
         i+=1;
      }
      assert_eq!( i,9, "A total of 54 fragments!")

   }

}