#ifndef MagicCarpetDefinitions_h
#define MagicCarpetDefinitions_h

#define MAXVOICES 32           // maximum number of voices



//#define MAXTABLELENGTH 262144    // =2^18 
#define MAXTABLELENGTH 524288    // =2^19
//#define MAXTABLELENGTH 1048576   // =2^20

//#define MAXTABLELENGTH 8192  // for test only
//#define MAXTABLELENGTH 2048  // for test only
//#define MAXTABLELENGTH 256  // for test only

 enum sampleSlots
 {
  TOP_LEFT = 0,
  TOP_RIGHT,
  BOTTOM_LEFT,
  BOTTOM_RIGHT,
  X_MOD,
  Y_MOD
 };

#endif // #ifndef MagicCarpetDefinitions_h