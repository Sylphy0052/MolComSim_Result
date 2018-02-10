# ActiveとPassiveの差について

## excelについて
平均値，中央値，TxRx平均値についての差をグラフにしたもの.

まず，ほとんどがACTIVEの方が早い.
### Mean
距離が短く，duplicationが小さい(1を除く)ほど，平均値の差は大きくなる...
気はするが，規則性はほとんどないと思われる

### Median
中央値に関しては20%以上差があるものが目立った．

duplicationが1のものはすべて0~5の範囲で収まっている．

距離が短いほど，差は大きくなっている．

duplicationが30または50の時，差が一番大きくなる．そこから離れるほど差は小さくなる．

差が一番大きいものよりduplicationが大きい時，duplicationに差の大きさが反比例している．

### TxRx-Mean
ほとんどの結果が10以下になり，初めてPASSIVEの方が早い時が出てきた．

この結果は距離が大きいほど，duplicationが大きいほど，PASSIVEの方が早い.
->duplicationと距離が大きいとmicrotubuleが上手く機能しない?(衝突しやすくなるため)
->なぜACTIVEがPASSIVEより遅くなるのかはわからない

## 以下の図について
左がActive,右がPassive
上の方が値が小さい(duplicationかdistanceか)

## Cumprobについて(compare_cumprob_each_duplication.png参照)
大した差はなかった

## Jitterについて(compare_jitter_by_duplication_each_distance.png参照)
こちらも大した差はなかった
duplicationが1->10で一番差が大きい
10->20も半分くらいになっている
そのあとはほとんど同じ

## Meanについて(compare_mean_by_diffusioncoefficient.png参照)
TxRxMeanとAnalyticalを比較すると，ACTIVEはTxRxMeanの方が大きく，PASSIVEはAnalyticalの方が大きい
-> Analyticalモデルの計算方法の違いによるもの?
Analyticalモデルは距離が小さいと精度が高い(距離30の時はほとんど同じ)

## Meanについて(compare_mean_by_distance_each_duplication.png参照)
グラフの形はほとんど同じ

## Meanについて(compare_mean_by_duplication_each_distance.png参照)
グラフの形はほとんど同じ

## Medianについて(compare_median_by_distance_each_duplication.png参照)
グラフの形はほとんど同じ

## TxRxMeanについて(compare_txrx_mean_by_distance_each_duplication.png参照)
PASSIVEは距離が大きくなっても近い値になっているが，ACTIVEは距離が大きくなるとばらつきが大きくなっている
ACTIVEはAnalytical RTTと大きく離れている -> 計算式の問題?

## jitterの回帰式について(regression_jitter_by_duplication_each_distance.png参照)
グラフの形も，方程式の値もほとんど同じ

## Medianの回帰式について(regression_median_by_duplication_each_distance.png参照)
グラフの形も，方程式の値もほとんど同じ

そもそもjitterの回帰式とMedianの回帰式もほとんど同じ
