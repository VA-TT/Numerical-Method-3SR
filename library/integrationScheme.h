#ifndef INTEGRATION_SCHEME_H
#define INTEGRATION_SCHEME_H

#include <DualDifferentiation.h>

template <typename T> void euler(T &pos, T &vel, const T &acc, const T &dt) {
  pos += dt * vel;
  vel += dt * acc;
}

template <typename T> void modiEuler(T &pos, T &vel, const T &acc, T dt) {
  // current acceleration
  //   T acc = F(pos,vel)/m;

  // Predictor
  T pos_predict = pos + vel * dt;
  T vel_predict = vel + acc * dt;

  // recompute acceleration at predicted state
  //   T acc_predict = F(pos_predict, vel_predict)/m;

  // Corrector
  pos += dt / T{2} * (vel + vel_predict);
  vel += dt / T{2} * (acc + acc_predict);
}

#endif