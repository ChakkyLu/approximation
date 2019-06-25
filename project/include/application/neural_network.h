#ifndef _APPLICATION_NN_H_INCLUDED
#define _APPLICATION_NN_H_INCLUDED

#include "global.h"

namespace neuralnetwork {

  typedef std::vector<float> singleVector;
  typedef std::vector<float> weightVector;
  // typedef std::vector<std::vector<int>> doubleVectore;
  typedef std::vector<std::vector<float>> doubleVectore;
  typedef std::vector<int> labelVector;

  class Network {
  public:
    Network() {};

    virtual ~Network() { Clear(); };

    void Clear() {
    }

  public:
    void init(float bias=1, int epoch=30, float lr=0.1, bool early_stopping=false);
    void fit(doubleVectore inputs, labelVector label);
    // void fit(doubleVectore inputs, singleVector label);
    float predict(doubleVectore inputs);
    int evaluate(doubleVectore test_data);
    int dotMultiply(singleVector v1, weightVector v2);
    int predict(singleVector inputs);


  public:
    weightVector weights;
    float bias;
    int epoch;
    float lr;
    bool early_stopping;
  };
}

#endif _APPLICATION_NN_H_INCLUDED
