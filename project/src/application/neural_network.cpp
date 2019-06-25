#include "neural_network.h"
#include <vector>

using namespace std;

namespace neuralnetwork {

  int Network::dotMultiply(singleVector v1, weightVector v2) {
    long int size = v1.size();
    if (size != v2.size()) {
      cout << "Matrix size does not match, cannot compute" << endl;
      return 0;
    }
    float product = 0;
    for (long int i=0; i<size; i++) {
      product += v1[i] * v2[i];
    }
    return (int)product;
  }

  void Network::init(float bias, int epoch, float lr, bool early_stopping) {
    Network::bias = bias;
    Network::epoch = epoch;
    Network::lr = lr;
    Network::early_stopping =early_stopping;
  }

  int Network::predict(singleVector inputs) {
    inputs.push_back(0.0);
    return Network::dotMultiply(inputs, Network::weights);
  }

  void Network::fit(doubleVectore inputs, labelVector label) {
    Network::weights.resize(inputs[0].size()+1, 0.0);
    singleVector tmp_vector;
    int real_y = 0;
    int test_y = 0;

    cout << "Trainning the Model" << endl;

    for (int i=0; i<Network::epoch; i++) {
      double cost = 0;
      cout << "Epoch " << i << "\t";
      for (long int j=0; j<inputs.size(); j++) {
        tmp_vector.clear();
        tmp_vector = inputs[j];
        tmp_vector.push_back(Network::bias);
        real_y = label[j];
        test_y = Network::dotMultiply(tmp_vector, Network::weights);
        if (test_y != real_y) {
          int flag = test_y < real_y ? 1 : -1;
          for (long int k=0; k<Network::weights.size(); k++) {
            Network::weights[k] += flag * tmp_vector[k] * Network::lr;
          }
        }
      }
      for (long int j=0; j<inputs.size(); j++) {
        if (Network::predict(inputs[j]) != label[j]) {
          cost += 1;
        }
      }
      cout << "Cost: " << cost / inputs.size() << endl;
    }
  }

}
