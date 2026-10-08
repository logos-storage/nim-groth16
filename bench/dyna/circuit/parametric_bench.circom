
// 
// parametric circuit for benchmarking dynamic proofs
// 

pragma circom 2.1.1;

//------------------------------------------------------------------------------
// The S-box

template SBox() {
  signal input  inp;
  signal output out;

  signal x2 <== inp*inp;
  signal x4 <== x2*x2;

  out <== inp*x4;
}

//------------------------------------------------------------------------------
// full round

template ExternalRound(i) {
  signal input  inp[3];
  signal output out[3];

  var round_consts[8][3] =

    [ [ 0x2c4c51fd1bb9567c27e99f5712b49e0574178b41b6f0a476cddc41d242cf2b43
      , 0x1c5f8d18acb9c61ec6fcbfcda5356f1b3fdee7dc22c99a5b73a2750e5b054104
      , 0x2d3c1988b4541e4c045595b8d574e98a7c2820314a82e67a4e380f1c4541ba90
      ]
    , [ 0x052547dc9e6d936cab6680372f1734c39f490d0cb970e2077c82f7e4172943d3
      , 0x29d967f4002adcbb5a6037d644d36db91f591b088f69d9b4257694f5f9456bc2
      , 0x0350084b8305b91c426c25aeeecafc83fc5feec44b9636cb3b17d2121ec5b88a
      ]
    , [ 0x1815d1e52a8196127530cc1e79f07a0ccd815fb5d94d070631f89f6c724d4cbe
      , 0x17b5ba882530af5d70466e2b434b0ccb15b7a8c0138d64455281e7724a066272
      , 0x1c859b60226b443767b73cd1b08823620de310bc49ea48662626014cea449aee
      ]
    , [ 0x1b26e7f0ac7dd8b64c2f7a1904c958bb48d2635478a90d926f5ff2364effab37
      , 0x2da7f36850e6c377bdcdd380efd9e7c419555d3062b0997952dfbe5c54b1a22e
      , 0x17803c56450e74bc6c7ff97275390c017f682db11f3f4ca6e1f714efdfb9bd66
      ]
    , [ 0x25672a14b5d085e31a30a7e1d5675ebfab034fb04dc2ec5e544887523f98dede
      , 0x0cf702434b891e1b2f1d71883506d68cdb1be36fa125674a3019647b3a98accd
      , 0x1837e75235ff5d112a5eddf7a4939448748339e7b5f2de683cf0c0ae98bdfbb3
      ]
    , [ 0x1cd8a14cff3a61f04197a083c6485581a7d836941f6832704837a24b2d15613a
      , 0x266f6d85be0cef2ece525ba6a54b647ff789785069882772e6cac8131eecc1e4
      , 0x0538fde2183c3f5833ecd9e07edf30fe977d28dd6f246d7960889d9928b506b3
      ]
    , [ 0x07a0693ff41476abb4664f3442596aa8399fdccf245d65882fce9a37c268aa04
      , 0x11eb49b07d33de2bd60ea68e7f652beda15644ed7855ee5a45763b576d216e8e
      , 0x08f8887da6ce51a8c06041f64e22697895f34bacb8c0a39ec12bf597f7c67cfc
      ]
    , [ 0x2a912ec610191eb7662f86a52cc64c0122bd5ba762e1db8da79b5949fdd38092
      , 0x2031d7fd91b80857aa1fef64e23cfad9a9ba8fe8c8d09de92b1edb592a44c290
      , 0x0f81ebce43c47711751fa64d6c007221016d485641c28c507d04fd3dc7fba1d2
      ]
    ];

  component sb[3];
  for(var j=0; j<3; j++) {
    sb[j] = SBox();
    sb[j].inp <== inp[j] + round_consts[i%8][j];
  }

  out[0] <== 2*sb[0].out +   sb[1].out +   sb[2].out;
  out[1] <==   sb[0].out + 2*sb[1].out +   sb[2].out;
  out[2] <==   sb[0].out +   sb[1].out + 3*sb[2].out;
}

//------------------------------------------------------------------------------
// Parametric permutation

template Permutation(nrounds) {
  signal input  inp[3];
  signal output out[3];

  signal aux[nrounds+1][3];

  aux[0] <== inp;

  component ext[nrounds];
  for(var k=0; k<nrounds; k++) { ext[k] = ExternalRound(k); }

  for(var k=0; k<nrounds; k++) {
    ext[k].inp <== aux[k  ];
    ext[k].out ==> aux[k+1];
  }

  out <== aux[nrounds];
}

//------------------------------------------------------------------------------
// Parametric hash

template Hash(nrounds) {
  signal input  inp;
  signal output out;

  component perm = Permutation(nrounds);

  perm.inp[0] <== inp;
  perm.inp[1] <== inp+1;
  perm.inp[2] <== inp+2;

  out <== perm.out[0] + perm.out[1] + perm.out[2];
}

//------------------------------------------------------------------------------
// Parametric benchmark
//
// expected number of constraints = 4 + 9*( staticRounds + updateRounds )
// 
// where the 4 is the public IO malleability equations 
// (not included in the circom-compiled R1CS)
//

template Benchmark(staticRounds,updateRounds) {

  signal input  staticInp;
  signal input  updateInp;
  signal input  secretInp;
  signal output out;

  component staticComp = Hash(staticRounds);
  component updateComp = Hash(updateRounds);

  staticComp.inp <== staticInp + secretInp;
  updateComp.inp <== updateInp + staticComp.out;
  updateComp.out ==> out;
} 

//------------------------------------------------------------------------------
