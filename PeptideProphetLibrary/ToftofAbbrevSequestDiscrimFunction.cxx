#include "ToftofAbbrevSequestDiscrimFunction.h"

ToftofAbbrevSequestDiscrimFunction::ToftofAbbrevSequestDiscrimFunction() : SequestDiscrimFunction(0) {

  float consts[] = {-1.244};
  float xcorrs[] = {5.375};
  float deltas[] = {5.078};
  float ranks[] = {-0.166};
  float massdiffs[] =  {0.0};

  const_ = consts[charge_];
  xcorr_p_wt_ = xcorrs[charge_];
  delta_wt_ = deltas[charge_];
  log_rank_wt_ = ranks[charge_];
  abs_massd_wt_ = massdiffs[charge_];

}
