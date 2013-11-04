#include <boost/shared_ptr>

class IConfigReader
{
  ReadTomoParameters(float_tt &tomoStart, float_tt &tomoStep, int &tomoCount,
                     float_tt &zoomFactor);
  ReadAberrationAmplitudes(float_tt &Cs, float_tt &C5,
                           float_tt &df0, float_tt &astig,
                           float_tt &a33, float_tt &a31,
                           float_tt &a44, float_tt &a42,
                           float_tt &a55, float_tt &a53, float_tt &a51,
                           float_tt &a66, float_tt &a64, float_tt &a62);
  ReadAberrationAngles(float_tt &astig,
                       float_tt &phi33, float_tt &phi31,
                       float_tt &phi44, float_tt &phi42,
                       float_tt &phi55, float_tt &phi53, float_tt &phi51,
                       float_tt &phi66, float_tt &phi64, float_tt &phi62);
};

typedef boost::shared_ptr<IConfigReader> ConfigReaderPtr






