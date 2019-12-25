//#include "rosic_BeatDetector.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

BeatDetector::BeatDetector()
{
  
}

BeatDetector::~BeatDetector()
{
  
}


//-------------------------------------------------------------------------------------------------
// setup:




//-------------------------------------------------------------------------------------------------
// inqiury:




//-------------------------------------------------------------------------------------------------
// processing:

void BeatDetector::processOnsets()
{
  estimateTempi();
  computeBeatProbabilities();
  markBeats();
}

void BeatDetector::estimateTempi()
{
  if( onsets.size() == 0 )
    return;

  // algorithm parameters:
  float              minBPM   = 90.f;     // minimum expected tempo in BPM
  float              margin   = 0.125f;   // relative tolerance margin for the clusters 
  unsigned int       M        = 20;       // memory (number of IOIs to consider)
  const unsigned int numTempi = 5;        // number of potential tempi to maintain - use size_t

  // internal variables:
  float        maxBeatInterval = 60.f/minBPM;
  float        minBeatInterval = maxBeatInterval/2.f;
  unsigned int numOnsets       = (unsigned int) onsets.size(); // use size_t

  // allocate memory for the tempo and confidence estimates:
  float** tempoEstimates = new float*[numOnsets];
  float** confidences    = new float*[numOnsets];
  for(unsigned int o=0; o<numOnsets; o++)
  {
    tempoEstimates[o] = new float[numTempi];
    confidences[o]    = new float[numTempi];
  }

  for(unsigned int o=0; o<numOnsets; o++)
  {
    for(unsigned int t=0; t<numTempi; t++)
    {
      tempoEstimates[o][t] = 0.f;
      confidences[o][t]    = 0.f;
    }
  }

  // cluster the onsets:
  std::vector<float> clusterMeans;  
  std::vector<float> clusterWeights;
  for(unsigned int n=1; n<numOnsets; n++)
  {
    unsigned int m = RAPT::rsMin(n, M);   // number of IOIs to consider

    // cluster the onsets:
    clusterMeans.clear();
    clusterWeights.clear();

    for(unsigned int k=0; k<m; k++)
    {
      float ioi = (onsets[n-k].timeInSamples - onsets[n-k-1].timeInSamples) / (float) sampleRate;
      ioi       = wrapIntoOctave(ioi, minBeatInterval);
      float w   = 1.f;  // ioi weight   

      // alternative formulas for the ioi-weight:
      //float w = onsets[n-k].strength * onsets[n-k-1].strength;
      //float w = onsets[n-k].strength + onsets[n-k-1].strength;

      // check if the current ioi is within the tolerance margin of some existing cluster - if so, 
      // add it to the cluster and update the cluster's mean:  
      bool clusterFound = false; 
      for(unsigned int c=0; c<clusterMeans.size(); c++)
      {
        if( fabs(clusterMeans[c]-ioi) < margin*clusterMeans[c] )
        {
          // update the weighted running mean of this cluster:
          float cw     = clusterWeights[c];          // old cluster weight
          float rwm    = clusterMeans[c];            // old running cluster mean    
          float cwNew  = clusterWeights[c] + w;      // new cluster weight
          float rwmNew = (cw*rwm + w*ioi) / cwNew;   // new running cluster mean  

          clusterMeans[c]   = rwmNew;
          clusterWeights[c] = cwNew;

          clusterFound = true;
          break;
        }
      }
      if( clusterFound == false || clusterMeans.size() == 0 )
      {
        // create new cluster with mean equal to the new ioi:
        clusterMeans.push_back(ioi);
        clusterWeights.push_back(w);
      }
    } // end of clustering

    // find the 'numTempi' top-rated tempo estimates and write them into the tempo-tracks:
    float weightSum = 0.0;
    for(unsigned int c=0; c<clusterWeights.size(); c++)
      weightSum += clusterWeights[c];

    for(unsigned int t=0; t<RAPT::rsMin(numTempi,(unsigned int) clusterWeights.size()); t++)
    {
      // find cluster with maximum weight:
      float maxWeight = 0.f;
      int   maxIndex  = 0;
      for(unsigned int i=0; i<clusterWeights.size(); i++)
      {
        if( clusterWeights[i] > maxWeight )
        {
          maxIndex  = i;
          maxWeight = clusterWeights[i];
        }
      }

      // write the cluster's data into the matrix of tempo-estimates (and associated confidences):
      tempoEstimates[n][t]     = 60.f/clusterMeans[maxIndex];
      confidences[n][t]        = maxWeight/weightSum;
      clusterWeights[maxIndex] = 0.f;  // to not find the same maximum weight again
    }
  }

  //---------------------------------------------------------------------------
  // reorder tempo-columns:

  // obtain the median of each column in the tempo matrix (each column of this matrix represents a 
  // time varying tempo estimate):
  float medianTempi[numTempi];
  float *tempoTrack = new float[numOnsets];
  for(unsigned int t=0; t<numTempi; t++)
  {
    // find median:
    for(unsigned int o=0; o<numOnsets; o++)
      tempoTrack[o] = tempoEstimates[o][t];
    medianTempi[t] = RAPT::rsArrayTools::median(tempoTrack, numOnsets);
  }
  delete[] tempoTrack;

  //float tempoEstimates2[16][5];   // \todo: allocate this memory with new
  //float confidences2[16][5];      // \todo: allocate this memory with new
  float** tempoEstimates2 = new float*[numOnsets];
  float** confidences2    = new float*[numOnsets];
  for(unsigned int o=0; o<numOnsets; o++)
  {
    tempoEstimates2[o] = new float[numTempi];
    confidences2[o]    = new float[numTempi];
  }

  // re-order the values in the matrix so as to obtain consistency along the colums - in each row 
  // we select the value from the column that is closest to the respective column-median
  for(unsigned int o=0; o<numOnsets; o++)
  {
    for(unsigned int m=0; m<numTempi; m++)
    {
      float med     = medianTempi[m];
      int   index   = 0;
      float minDiff = FLT_MAX;

      for(unsigned int t=0; t<numTempi; t++)
      {
        if( fabs(tempoEstimates[o][t]-med) < minDiff )
        {
          index   = t;
          minDiff = fabs(tempoEstimates[o][t]-med);
        }
      }
      tempoEstimates2[o][m] = tempoEstimates[o][index];
      confidences2   [o][m] = confidences   [o][index];
    }
  }

  // the first colums now contains hopefully the right tempo-estimate for each onset - copy the 
  // values into the onset's bpm members:
  onsets[0].bpm = tempoEstimates2[1][0]; // first onset must be treated seperately by copying the 
                                        // value from the second one 
  for(unsigned int o=1; o<numOnsets; o++)
    onsets[o].bpm = tempoEstimates2[o][0];

  // clean up temporarily allocated memory:
  for(unsigned int o=0; o<numOnsets; o++)
  {
    delete[] tempoEstimates[o];
    delete[] confidences[o];
    delete[] tempoEstimates2[o];
    delete[] confidences2[o];
  }
  delete[] tempoEstimates;
  delete[] confidences;
  delete[] tempoEstimates2;
  delete[] confidences2;
}

void BeatDetector::computeBeatProbabilities()
{
  if( onsets.size() == 0 )
    return;
  int numOnsets = (int) onsets.size();

  // algorithm parameters:
  int   searchSpan = 8;    // number of beat-intervals, over which we sum the coincidences of onsets 
                           // and beats
  int   numPasses  = 8;    // number of passes of the whole loop
  float tolerance  = 0.1f; // maximum relative time deviation for coincidence

  float *beatProbabilitiesTmp = new float[numOnsets];
  float *beatProbabilities    = new float[numOnsets];
  for(int o=0; o<numOnsets; o++)
  {
    beatProbabilities[o]    = onsets[o].strength;
    beatProbabilitiesTmp[o] = onsets[o].strength;
  }

  // outer loop over the number of passes
  for(int pass=1; pass<=numPasses; pass++)
  {
    RAPT::rsArrayTools::copy(beatProbabilities, beatProbabilitiesTmp, numOnsets);
    for(int o=0; o<numOnsets; o++)
    {
      int   t  = onsets[o].timeInSamples;         // onset time in samples
      float bi = 60.f*sampleRate / onsets[o].bpm; // beat interval in samples
      float normalizer = 1.f;                     // to take into account boundary problems

      for(int j=1; j<=searchSpan; j++)
      {
        float tp = t + j*bi;  // forward predicted beat/onset time
        float to;             // for the observed onset time

        // search for an actual onset near the predicted onset:
        int k = o+j;
        while( true )
        {
          if( k >= numOnsets )
          {
            k = -1;   // indicates, that no onset was found
            break;
          }

          to = (float) onsets[k].timeInSamples; // observed onset time at index k

          if( fabs(tp-to) < tolerance*bi )
            break;    // nearby onset found at index k - leave while loop
          else if( (to-tp) >= tolerance*bi ) // no abs here
          {
            k = -1;   // we are already too far - indicate, that no onset was found  
            break;
          }
          k          = k+1;
          normalizer = normalizer+1;
        }
  
        // if some onset was found near the predicted beat, increase the beat-probability of the 
        // onset under investigation by an amount equal to the (preliminary) beat-probability of 
        // the found onset:
        if( k != -1 )
          beatProbabilities[o] = beatProbabilities[o] + beatProbabilitiesTmp[k];

        //...and now do the same procedure with backward prediction:
        tp = t-j*bi;
        k  = o-j;
        while( true )
        {
          if( k < 0 )
          {
            k = -1;
            break;
          }
          to = (float) onsets[k].timeInSamples;
          if( fabs(tp-to) < tolerance*bi )
            break;
          else if( (tp-to) >= tolerance*bi )
          {
            k = -1;
            break;
          }
          k          = k-1;
          normalizer = normalizer+1;    
        }
        if( k != -1 )
          beatProbabilities[o] = beatProbabilities[o] + beatProbabilitiesTmp[k];
      }

      beatProbabilities[o] /= normalizer;
    }
  }

  for(int o=0; o<numOnsets; o++)
    onsets[o].beatProbability = beatProbabilities[o];

  delete[] beatProbabilitiesTmp;
  delete[] beatProbabilities;
}

void BeatDetector::markBeats()
{
  if( onsets.size() == 0 )
    return;
  int numOnsets = (int) onsets.size();

  // algorithm parameters:
  float tolerance = 0.25;

  // set all the beat flags preliminarily to true - some will be falsified later:
  for(int o=0; o<numOnsets; o++)
    onsets[o].isBeat = true;

  for(int o=0; o<numOnsets; o++)
  {
    float t  = (float) onsets[o].timeInSamples;           // onset time in samples
    float bi = (float) (60.0*sampleRate / onsets[o].bpm); // beat interval in samples 

    // set the beat flag of onset i to zero, if an onset with higher beat-probability is found 
    // within a range slightly less than a beat-interval surrounding the onset o:
    float d = (1-tolerance) * bi;
    int   j = 0;
    while( (o+j) < numOnsets && onsets[o+j].timeInSamples-t < d )
    {
      if( onsets[o+j].beatProbability > onsets[o].beatProbability )
        onsets[o].isBeat = false;
      j++;
    }
      
    j = 0;
    while( (o-j) >= 0 && t-onsets[o-j].timeInSamples < d )
    {
      if( onsets[o-j].beatProbability > onsets[o].beatProbability )
        onsets[o].isBeat = false;
      j++;
    }
  }
}

float BeatDetector::wrapIntoOctave(float x, float xMin)
{
  float y = x;
  while( y < xMin )
    y *= 2.f;
  while( y >= 2.f*xMin )
    y *= 0.5f;

  return y;
}






