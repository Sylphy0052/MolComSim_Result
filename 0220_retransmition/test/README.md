# MolComSim_Analyzer
## How to use
1. Make "dat" directory here.
2. Put config files in "dat" directory.
3. Make "result" directory here.
4. Put result datas in "result" directory.
5. Run "./execute_analysis.sh".
6. If you want to remove files made by "execute_analysis.sh", Run "./remove_result.sh".

## What are made by this program
1. index.html: 設定ファイル，出力値，シミュレーションの統計が見れる．
2. compare_cumprob_each_duplication: 重複度ごとのCumulativeProbability
3. compare_jitter_by_duplication_each_distance: 距離ごとの標準偏差
4. compare_mean_by_distance: 平均，予測値，TxRx平均
5. compare_mean_by_distance_each_duplication: 重複度ごとの平均
6. compare_mean_by_duplication_each_distance: 距離ごとの平均
7. compare_median_by_distance_each_duplication: 重複度ごとの中央値
8. compare_txrx_mean_by_distance_each_duplication: 重複度ごとのTxRx平均
9. regression_jitter_by_duplication_each_distance: 距離ごとの標準偏差の回帰
10. regression_med_by_duplication_each_distance: 距離ごとの中央値の回帰
