//
// Created by andre on 19/07/2022.
//

#ifndef SHA_3_UTILS_H
#define SHA_3_UTILS_H

/*
 * get index X in convention
 * */
std::vector<size_t> transformX{2,3,4,0,1};
/*
 * get index Y in convention
 * */
std::vector<size_t> transformY{2,1,0,4,3};

/* SHA-3 constants */
constexpr uint64_t RC[] = {
    0x0000000000000001L, 0x0000000000008082L, 0x800000000000808aL,
    0x8000000080008000L, 0x000000000000808bL, 0x0000000080000001L,
    0x8000000080008081L, 0x8000000000008009L, 0x000000000000008aL,
    0x0000000000000088L, 0x0000000080008009L, 0x000000008000000aL,
    0x000000008000808bL, 0x800000000000008bL, 0x8000000000008089L,
    0x8000000000008003L, 0x8000000000008002L, 0x8000000000000080L,
    0x000000000000800aL, 0x800000008000000aL, 0x8000000080008081L,
    0x8000000000008080L, 0x0000000080000001L, 0x8000000080008008L
};

uint64_t ROTL(const uint64_t &x, const size_t &n) {
  return (x << n) | (x >> (64 - n));
}

uint64_t ROTR(const uint64_t &x, const size_t &n) {
  return (x >> n) | (x << (64 - n));
}

std::vector<std::vector<u_int64_t>> theta(const std::vector<std::vector<u_int64_t>> &state){
  std::vector<std::vector<u_int64_t>> statePrime(5, std::vector<u_int64_t>(5, 0));
  std::vector<u_int64_t> C(5, 0);
  std::vector<u_int64_t> D(5, 0);

  for(size_t indexOfX = 0; indexOfX < 5; ++indexOfX){
    C[transformX[indexOfX]] = state[transformX[indexOfX]][0] ^ state[transformX[indexOfX]][1] ^ state[transformX[indexOfX]][2] ^ state[transformX[indexOfX]][3] ^ state[transformX[indexOfX]][4];
  }

  /* How is D calculated? (book convention)
   * C[x-1 mod 5] = {0, 1,2,3,4,5,6,7,8,9,10,...,61,62,63}
   * Xored with:     |  | | | | | | | | |  |     |  |  |
   * C[x+1 mod 5] = {63,0,1,2,3,4,5,6,7,8, 9,...,60,61,62}
   * Esto quiere decir que C[x+1 mod 5] tuvo una rotación circular a la derecha por 1 posición
   * */
  for(size_t indexOfX = 0; indexOfX < 5; ++indexOfX){
    D[transformX[indexOfX]] = C[transformX[(indexOfX - 1)%5]] ^ ROTR(C[transformX[(indexOfX + 1)%5]], 1);
  }

  for(size_t indexOfX = 0; indexOfX < 5; ++indexOfX){
    for(size_t indexOfY = 0; indexOfY < 5; ++indexOfY){
      statePrime[transformX[indexOfX]][transformY[indexOfY]] = state[transformX[indexOfX]][transformY[indexOfY]] ^ D[transformX[indexOfX]];
    }
  }

  return statePrime;
}

std::vector<std::vector<u_int64_t>> rho(const std::vector<std::vector<u_int64_t>> &state){
  std::vector<std::vector<u_int64_t>> statePrime(5, std::vector<u_int64_t>(5, 0));
  std::vector<u_int64_t> rotationsForT {1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 2, 14, 27, 41, 56, 8, 25, 43, 62, 18, 39, 61, 20, 44};

  statePrime[transformX[0]][transformY[0]] = state[transformX[0]][transformY[0]];

  size_t x = 1;
  size_t y = 0;
  for (size_t t = 0; t < 24; ++t) {
    statePrime[transformX[x]][transformY[y]] = ROTR(state[transformX[x]][transformY[y]], rotationsForT[t]);
    size_t tmp = x;
    x = y;
    y = (2*tmp+3*y)%5;
  }

  return statePrime;
}

std::vector<std::vector<u_int64_t>> pi(const std::vector<std::vector<u_int64_t>> &state){
  std::vector<std::vector<u_int64_t>> statePrime(5, std::vector<u_int64_t>(5, 0));

  for(size_t indexOfX = 0; indexOfX < 5; ++indexOfX){
    for(size_t indexOfY = 0; indexOfY < 5; ++indexOfY){
      statePrime[transformX[indexOfX]][transformY[indexOfY]] = state[transformX[(indexOfX + 3*indexOfY)%5]][transformY[indexOfX]];
    }
  }

  return statePrime;
}

std::vector<std::vector<u_int64_t>> chi(const std::vector<std::vector<u_int64_t>> &state){
  std::vector<std::vector<u_int64_t>> statePrime(5, std::vector<u_int64_t>(5, 0));

  for(size_t indexOfX = 0; indexOfX < 5; ++indexOfX){
    for(size_t indexOfY = 0; indexOfY < 5; ++indexOfY){
      statePrime[transformX[indexOfX]][transformY[indexOfY]] = state[transformX[indexOfX]][transformY[indexOfY]] ^ ((~state[transformX[(indexOfX + 1)%5]][indexOfY] ^ 1) & state[transformX[(indexOfX  + 2)%5]][transformY[indexOfY]]);
    }
  }

  return statePrime;
}

std::vector<std::vector<u_int64_t>> iota(const std::vector<std::vector<u_int64_t>> &state, const size_t i){
  std::vector<std::vector<u_int64_t>> statePrime(state);

  statePrime[transformX[0]][transformY[0]] ^= RC[i];

  return statePrime;
}

std::vector<std::vector<u_int64_t>> Rnd(const std::vector<std::vector<u_int64_t>> &state, const size_t i){

  return iota(chi(pi(rho(theta(state)))), i);
}

std::string KECCAK(const std::string &S){
  if(S.length() != 200)
    throw std::runtime_error("La cadena no tiene longitud de 1600 bits");

  std::string SPrime;
  std::vector<std::vector<u_int64_t>> state(5, std::vector<u_int64_t>(5, 0));

  /* Assign string to state array 3.1.2 */
  size_t z = 0;
  for(size_t indexOfY = 0; indexOfY < 5; ++indexOfY){
    for(size_t indexOfX = 0; indexOfX < 5; ++indexOfX){
      /* Fill 64 bit with 8 bits */
      uint64_t S_1 = S[z];
      uint64_t S_2 = S[z + 1];
      uint64_t S_3 = S[z + 2];
      uint64_t S_4 = S[z + 3];
      uint64_t S_5 = S[z + 4];
      uint64_t S_6 = S[z + 5];
      uint64_t S_7 = S[z + 6];
      uint64_t S_8 = S[z + 7];

      uint64_t S_z = (S_1 << 56) | (S_2 << 48) | (S_3 << 40) | (S_4 << 32) | (S_5 << 24) | (S_6 << 16) | (S_7 << 8) | S_8;

      /* String to state array */
      // Turn it to little endian
      state[transformX[indexOfX]][transformY[indexOfY]] = __builtin_bswap64(S_z);

      z += 8;
    }
  }

  /* Round */
  for (int i = 0; i < 24; ++i) {
    state = Rnd(state, i);
  }

  /* Convert state back to string */
  for(size_t indexOfY = 0; indexOfY < 5; ++indexOfY){
    for(size_t indexOfX = 0; indexOfX < 5; ++indexOfX){
      /* Convert 64 bits to uint8_t */
      uint64_t lane = state[transformX[indexOfX]][transformY[indexOfY]];
      uint8_t char_1 = lane >> 56;
      uint8_t char_2 = (lane << 8) >> 56;
      uint8_t char_3 = (lane << 16) >> 56;
      uint8_t char_4 = (lane << 24) >> 56;
      uint8_t char_5 = (lane << 32) >> 56;
      uint8_t char_6 = (lane << 40) >> 56;
      uint8_t char_7 = (lane << 48) >> 56;
      uint8_t char_8 = (lane << 56) >> 56;

      SPrime.push_back(char_1);
      SPrime.push_back(char_2);
      SPrime.push_back(char_3);
      SPrime.push_back(char_4);
      SPrime.push_back(char_5);
      SPrime.push_back(char_6);
      SPrime.push_back(char_7);
      SPrime.push_back(char_8);
    }
  }

  return SPrime;
}

std::string SPONGE(const std::string &N){
  const int b_in_bytes = 200; // 1600 bits ancho de f
  const int d_in_bytes = 64; // 512 bits
  const int r_in_bytes = 72; // 576 bits rate
  const int c_in_bytes = b_in_bytes - r_in_bytes; // 1024 bits capacidad

  std::string P(N.begin(), N.end());

  /* Padding */
  // 00000000
  char zero = 0x0;
  // 10000000
  char oneTrailingZeroes = (char) 0x80;
  // 00000001
  char trailingZeroesOne = 0x01;
  // 10000001
  char onetrailingZeroesOne = (char) 0x81;

  // Check if only remain one byte to multiple of 72
  if(P.length() % r_in_bytes == 71)
    P += onetrailingZeroesOne;
  // If remain more than one byte to multiple of 72
  else {

    P += oneTrailingZeroes;

    while((P.length()%r_in_bytes) < 71){
      P += zero;
    }

    P += trailingZeroesOne;
  }

  /* S and 0^c */
  std::string S(b_in_bytes, (char) 0);
  std::string cZeroes(c_in_bytes, (char) 0);

  /* Keccak */
  for(size_t indexOfP = 0; indexOfP < P.size(); indexOfP += 72){
    std::string completeBlock = P.substr(indexOfP, 72);
    completeBlock.insert(completeBlock.end(), cZeroes.begin(), cZeroes.end());
    std::string xoredArgument(completeBlock.size(), (char) 0);

    if(completeBlock.size() != b_in_bytes && S.size() != b_in_bytes && xoredArgument.size() != b_in_bytes)
      throw std::runtime_error("algo salió mal!! S o P completo no son el número de bytes deseado");

    for(size_t indexOfS = 0; indexOfS < S.size(); ++indexOfS){
      xoredArgument[indexOfS] = S[indexOfS] ^ completeBlock[indexOfS];
    }
    S = KECCAK(xoredArgument);
  }

  /* Get output */
  std::string Z;

  while(true){
    std::string Stmp(S.begin(), S.begin() + r_in_bytes);
    Z += Stmp;
    if(d_in_bytes <= Z.length()){
      break;
    } else {
      S = KECCAK(S);
    }
  }

  return Z.substr(0,d_in_bytes);
}

std::string SHA_3(const std::string &M){
  return SPONGE(M);
}

#endif //SHA_3_UTILS_H
