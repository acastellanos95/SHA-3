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
    D[transformX[indexOfX]] = C[(transformX[indexOfX] - 1)%5] ^ ROTR(C[(transformX[indexOfX] + 1)%5], 1);
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

  statePrime[transformX[0]][transformY[0]] = state[transformX[0]][transformY[0]];

  size_t x = transformX[1];
  size_t y = transformX[0];
  for (size_t t = 0; t < 24; ++t) {

  }

  return statePrime;
}

#endif //SHA_3_UTILS_H
