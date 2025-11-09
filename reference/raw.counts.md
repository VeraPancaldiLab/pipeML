# Raw counts

Raw gene expression matrix from bulk RNAseq.

## Usage

``` r
raw.counts
```

## Format

Matrix with genes as rows and samples as columns

## Source

Mariathasan et al. (2018), doi: https://doi.org/10.1038/nature25501

## Examples

``` r
data(raw.counts)
head(raw.counts)
#>         SAM7f0d9cc7f001 SAM4305ab968b90 SAMcf018fee2acd SAMcc4675f394a1
#> A1BG                 22               9              43            9821
#> NAT2                  0               0              13             173
#> ADA                 294             735             345             944
#> CDH2                115              24              57            4288
#> AKT3               2694             239            2194            1450
#> GAGE12F               0               0               0               0
#>         SAM49f9b2e57aa5 SAM2e7aa8fa0ab3 SAMdf3e42c8672a SAMd027124354ce
#> A1BG                 23              32              24              96
#> NAT2                  0               2               3               2
#> ADA                1123             492            1275            1562
#> CDH2                481             331             121            4106
#> AKT3                457             669            1376            3701
#> GAGE12F               0               0               0               0
#>         SAMe7bf6c015192 SAM6dd7ad1d797d SAM18039827e1b9 SAMc692536a795a
#> A1BG                 30             147              19              42
#> NAT2                  2               4               6               3
#> ADA                 994            2472             369             412
#> CDH2                443             966              59             285
#> AKT3               1342            2129            1765            2045
#> GAGE12F               0               0               0               0
#>         SAM9a2cf3c06fb3 SAM557dde1b9f3e SAM23aa15d4a0b0 SAM468a9e1dc821
#> A1BG                 22              31               9               9
#> NAT2                  7               0               1               0
#> ADA                1057             207             813            1331
#> CDH2                203             647             322             227
#> AKT3               4500            1931            1795             934
#> GAGE12F               0               0               0               0
#>         SAM81b71522417a SAMb963dda93cfd SAMbcbc7957c264 SAM7fb6987514a4
#> A1BG                 57              51              10              27
#> NAT2                  0               1               1               0
#> ADA                 913            2389             326            1299
#> CDH2               5402             565             496            2113
#> AKT3                876            1628            1398            2736
#> GAGE12F               0               0               0               0
#>         SAM63405b04ab2d SAM18bc1078bc15 SAMd1bd63734394 SAMe9ae8beb82fa
#> A1BG                 34              56              49              25
#> NAT2                  0               1               1              16
#> ADA                1795             281            2130             846
#> CDH2                879              65              85             758
#> AKT3               1308            1745             934            2286
#> GAGE12F               0               0               0               0
#>         SAMba7176afe070 SAMbe83eae4026e SAMe5bc41772bc9 SAM23095936e611
#> A1BG                 13              14            1923              18
#> NAT2                  0               0              65               0
#> ADA                 366             534            1361             721
#> CDH2                666             411            4627             282
#> AKT3                602             760            1147             709
#> GAGE12F               0               0               0               0
#>         SAM7114d99032ec SAMdb3f50c9129c SAMbf1a3ae828e6 SAM52e3fa3ad574
#> A1BG                 20              16              10              51
#> NAT2                  2               0              12               1
#> ADA                 512             232             718            1561
#> CDH2                240             297             241             247
#> AKT3                958            3776            4294            2567
#> GAGE12F               0               0               0               0
#>         SAMd4c0837b0997 SAM9cafb905b36a SAM032c642382a7 SAM0ce9c983b20f
#> A1BG                 11             144              38              22
#> NAT2                  2               4               1               4
#> ADA                 296            1187            2393            1493
#> CDH2                105            2884             181             202
#> AKT3               2031            2444            1471            1025
#> GAGE12F               0               0               0               0
#>         SAM2f228939632f SAM297c0301e861 SAM1fa6bcb7fc48 SAM6d2ae0c39b96
#> A1BG                 17              30              16               9
#> NAT2                 48               0               0               0
#> ADA                 435             621             447             926
#> CDH2                754           12878              62             149
#> AKT3               1227            5264             751             614
#> GAGE12F               0               0               0               0
#>         SAMeff2ce356ccb SAM110501d0eedb SAM2e9ac0b1b250 SAMc0ef41aa6c8b
#> A1BG                  8              37              94              82
#> NAT2                  0               0               1               0
#> ADA                 453             469             375             615
#> CDH2                 38              44            1268              41
#> AKT3                592             526            1632             602
#> GAGE12F               0               0               0               0
#>         SAMd636e3461955 SAM30b5c6c54cf7 SAMd98bac0a070f SAM8e8ef2368dfa
#> A1BG                 50              20              17              45
#> NAT2                  4               0               5               0
#> ADA                 880            1615             875             647
#> CDH2               1159             164             401            2310
#> AKT3               3139             443            1283            2707
#> GAGE12F               0               0               0               0
#>         SAMb4c7a001537d SAMfd947610629d SAM943df5cf15df SAM39eb94fa504d
#> A1BG                 11               7              15              91
#> NAT2                  1               2               1               2
#> ADA                 398             390            2797            1913
#> CDH2                 21             333             541             244
#> AKT3               2833            1380            1190             802
#> GAGE12F               0               0               0               0
#>         SAM62fb1388c871 SAMc0da5d48686d SAM2dc3f04e45e9 SAM4501e41e4751
#> A1BG                 22              33              17              29
#> NAT2                  4              10               6               4
#> ADA                 749            8599             956             864
#> CDH2                141             268             818            1090
#> AKT3               1711            6862            4964            5870
#> GAGE12F               0               0               0               0
#>         SAM5b57e47fdcb3 SAMae4da274eded SAMd35318127278 SAM75142fcab9df
#> A1BG                 13              30              19              50
#> NAT2                  1               3               0               3
#> ADA                 918            1212             988             980
#> CDH2                  6              90             127            1456
#> AKT3                772            1146             627            2772
#> GAGE12F               0               0               0               0
#>         SAM166a419a4e5a SAM025b45c27e05 SAMe07c4560772d SAM560f23d6a3ad
#> A1BG                 32              77               9              33
#> NAT2                  1               2               0               2
#> ADA                1011            2116             363            1397
#> CDH2                979            2451              12            1858
#> AKT3               4063            3509             535             592
#> GAGE12F               0               0               0               0
#>         SAMf3a9bce50099 SAMdab9ca8fb5de SAM34430ef08e5b SAM4b0175e8db6e
#> A1BG                 38              18              36              49
#> NAT2                  0               1               2               5
#> ADA                 514            1020             457             625
#> CDH2                190             414              48             493
#> AKT3               1121            3479            1461            1829
#> GAGE12F               0               0               0               0
#>         SAM54e58f1b0230 SAM5e3bae090b8c SAM3cb94b0d5297 SAMcb132b0cdd2c
#> A1BG                 10              18              18              32
#> NAT2                  4               0               0               0
#> ADA                  69             293             423            1027
#> CDH2                127             166             176             321
#> AKT3                707            1343            1557            5662
#> GAGE12F               0               0               0               0
#>         SAMe97af0feefdf SAMee3844cc0b9f SAMdad5c29dc105 SAMa424c75831b4
#> A1BG                 21              49              32              23
#> NAT2                  1               2               0               0
#> ADA                 516             836            1152             251
#> CDH2                192            1048             327             151
#> AKT3               1028            1967            1108            1467
#> GAGE12F               0               0               0               0
#>         SAM45c8e6412c66 SAM73663ee4a96e SAM0bdb3428bd13 SAM3330c03fdf00
#> A1BG                 20              30              22              34
#> NAT2                  1               1              16               4
#> ADA                3146             541             942            1868
#> CDH2               1317             185             197             588
#> AKT3               2021             771            1512            5115
#> GAGE12F               0               0               0               0
#>         SAM0a7c2091dd56 SAM5d1dfd5207f5 SAM99a46b9eec27 SAM1ac4e3dee297
#> A1BG                 74               7              57               9
#> NAT2                  0               2               7               3
#> ADA                3517             664            1183             214
#> CDH2                490              55            1802              24
#> AKT3                857            1221            5177             265
#> GAGE12F               0               0               0               0
#>         SAM553c3c35b847 SAMc2a1820d4e6b SAMabc151b01ea3 SAMe1eb5d988760
#> A1BG                 55              28            5325              17
#> NAT2                  4               1             135               7
#> ADA                1255             928             643             895
#> CDH2               1563             359            1766             305
#> AKT3               2290            2245             606            2677
#> GAGE12F               0               0               0               0
#>         SAMe3210d3632b4 SAMb15ad09d6e24 SAM7893196e0e89 SAM961d04c42bd9
#> A1BG                 24              32              24              17
#> NAT2                  0               1               0               3
#> ADA                1301            1162             887            1208
#> CDH2               1816             298             173             232
#> AKT3               2393            3368            2211            1074
#> GAGE12F               0               0               0               0
#>         SAMb0d11db9aa79 SAM19fec8f3b3bd SAMa90d73f8d891 SAM415f36ad349e
#> A1BG                 18              71              13              35
#> NAT2                  1               0               3               5
#> ADA                 590            1260            2510             689
#> CDH2                324             238             118             178
#> AKT3               1075            2007             728            2379
#> GAGE12F               0               0               0               0
#>         SAMad83c9c53537 SAM3ee5dcd894f0 SAM8a1b0e02ee42 SAMd7d57ee3a863
#> A1BG                 47              59             159              26
#> NAT2                  0               2              12               0
#> ADA                1622             487            1227             997
#> CDH2                547             394              52             750
#> AKT3               1414             592            2784            5669
#> GAGE12F               0               0               0               0
#>         SAMeb29625f76a5 SAM563d6233dfa2 SAMa1871f491b02 SAMd135d5867fe3
#> A1BG                 15              34               6              71
#> NAT2                  3               0              13               1
#> ADA                 492             189             417             537
#> CDH2                439             454             172             998
#> AKT3                906            1406             762            2789
#> GAGE12F               0               0               0               0
#>         SAM065890737112 SAMb470eb8f04be SAM675a12a09c15 SAM1e9c4d1d39ae
#> A1BG                 10              36              23              23
#> NAT2                  1               0               2               5
#> ADA                 368             396             863             793
#> CDH2                  9             134             606              36
#> AKT3                366            1494            1294             701
#> GAGE12F               0               0               0               0
#>         SAM4b7ea015fd9e SAM9306c5c92444 SAMd86389d0d768 SAM63b2189c36d7
#> A1BG                  4              64              19              11
#> NAT2                  4               3               2               0
#> ADA                 222             578              91             393
#> CDH2                 41              14              32              90
#> AKT3                354            4720             574             666
#> GAGE12F               0               0               0               0
#>         SAM18be5b395318 SAM0d855cff64e6 SAM6cbc10abddb0 SAMdee1011782cd
#> A1BG                  8              11               8              27
#> NAT2                  0               4               0               4
#> ADA                 209            1408             436            1031
#> CDH2                115              98             152             618
#> AKT3               1815             562             778            2155
#> GAGE12F               0               0               0               0
#>         SAM203dcf14f927 SAMe41b1e773582 SAM978a587b207e SAM5234688806a7
#> A1BG                 50               9              50               3
#> NAT2                  0               0               2               1
#> ADA                1322             459            1247             253
#> CDH2                120             161             200              43
#> AKT3               3856             844            2684             203
#> GAGE12F               0               0               0               0
#>         SAM2c9586161ce6 SAM76a431ba6ce1 SAMd3bd67996035 SAMd3601288319e
#> A1BG                  2              19              34              20
#> NAT2                  6               4               0               4
#> ADA                 143             458            1687             486
#> CDH2                 69             899             747             535
#> AKT3                499            1392            2547            4659
#> GAGE12F               0               0               0               0
#>         SAMba1a34b5a060 SAM18a4dabbc557 SAMfed609955db9 SAMf2aae1443f67
#> A1BG                 35              15              10              40
#> NAT2                  2               6               1               2
#> ADA                1699            1104             539            1368
#> CDH2                656              38             459             150
#> AKT3               1969             471            1617             978
#> GAGE12F               0               0               0               0
#>         SAM2bba8cb35e48 SAM5c139c5c1c4f SAM6780ed436b55 SAM85f0a3ac1c45
#> A1BG                 14              14              59               3
#> NAT2                  0               0               2               1
#> ADA                 195             423            2044             755
#> CDH2                562              32             362              64
#> AKT3               1862             702             716             359
#> GAGE12F               0               0               0               0
#>         SAM9d2494119c05 SAM5d989c86255e SAM49d48750e294 SAM7aa01fc49a80
#> A1BG               8097              20              22              35
#> NAT2                193               0               1               0
#> ADA                 772             465             832             938
#> CDH2               2289             772            1219             133
#> AKT3                418             920            2104            2791
#> GAGE12F               0               0               0               0
#>         SAMa321770ac31c SAM3894ac3956a5 SAM1c0ecfb3eb63 SAMbcb07ba81cee
#> A1BG                 26               9              27               3
#> NAT2                  0               3               0               0
#> ADA                1936             685            2692             357
#> CDH2                995             699             718               9
#> AKT3                508            2913            1822             374
#> GAGE12F               0               0               0               0
#>         SAM9aa6a095a9d6 SAM957378bd907f SAMf82bbdc267c8 SAM2624229effe8
#> A1BG                 31              46               7              23
#> NAT2                  1               0               0               2
#> ADA                1412            1277             147             595
#> CDH2                744             249              92             216
#> AKT3               1785            2730             925            1721
#> GAGE12F               0               0               0               0
#>         SAMd2492b2a31bb SAM670649e105b5 SAM5fe7a81a39dd SAM29da928587ad
#> A1BG                 53              24              12              32
#> NAT2                  7               0               2               0
#> ADA                1121             905             484             745
#> CDH2               1589             923             418             253
#> AKT3               1575            1527            1402             629
#> GAGE12F               0               0               0               0
#>         SAM0a0f2bac4b20 SAM8e469834acc1 SAM181b638b8248 SAM3b15b4c6311d
#> A1BG                 19              19              38              61
#> NAT2                  5               1               1               2
#> ADA                 898             279             675            3215
#> CDH2               1213             355            1648             988
#> AKT3               2947            2099            4676            3525
#> GAGE12F               0               0               0               0
#>         SAM1c8b086175ca SAMeaa477a5384b SAM1dda30f1c5be SAM9daccafc18db
#> A1BG                 25              23              51              21
#> NAT2                  0               0               0               0
#> ADA                 360             793             747             311
#> CDH2               1047              26            4983             684
#> AKT3               1872             539            5154            1083
#> GAGE12F               0               0               0               0
#>         SAM25510f300d79 SAM187e056d6a2a SAMbc8dc3a7b54e SAMae02629a97f7
#> A1BG                 41              27              32              21
#> NAT2                  1               1               0               1
#> ADA                 436             495             501             317
#> CDH2                619             192              42             275
#> AKT3                758             558             740             536
#> GAGE12F               0               0               0               0
#>         SAM7829a341b9f3 SAM0f956e757453 SAM09c84ec0cf34 SAMb2e4a082541a
#> A1BG                 25              20              13             152
#> NAT2                  1               2               0               0
#> ADA                 407             484             393            1433
#> CDH2                 81             155             408             457
#> AKT3                678            1647             739            2296
#> GAGE12F               0               0               0               0
#>         SAM771445e92421 SAM59fda9035d1d SAM6f2a102a99df SAMe56c96c51190
#> A1BG                 31              60               5              76
#> NAT2                  0               0               7               1
#> ADA                 225             555             427            1192
#> CDH2                431            3028             174            2411
#> AKT3               1500            1315            2422            2359
#> GAGE12F               0               0               0               0
#>         SAM1ab1b28d9f2b SAM6662f5181f87 SAM18b9351e265a SAM7fb7a13c096b
#> A1BG                 18              20              39              40
#> NAT2                  1               0              33               0
#> ADA                 138            1801            1673             272
#> CDH2               3877             526             117             305
#> AKT3                639             836             989             619
#> GAGE12F               0               0               0               0
#>         SAM87a8e18eb45b SAM9725303dce0c SAM6ff654a20f98 SAM94859b440b1d
#> A1BG                  2             146              22              13
#> NAT2                  0               3               1               0
#> ADA                1039            1153             435             276
#> CDH2                 15              41             210              65
#> AKT3                227             500             968             470
#> GAGE12F               0               0               0               0
#>         SAMc57eadb2d82b SAMaabf4afe4213 SAMce39dd79b441 SAM6964a6d7b967
#> A1BG                  8              35              17              19
#> NAT2                  0               0               1               4
#> ADA                  32            1661            1431             174
#> CDH2                 11             612             700              26
#> AKT3                538            2768             209             440
#> GAGE12F               0               0               0               0
#>         SAMf20b827dca51 SAM6792d6e98068 SAM31d9176e11fb SAM5cc2d9036053
#> A1BG                 37               9              18              24
#> NAT2                  2               0               0               0
#> ADA                1225              98             738             178
#> CDH2                 57              22             132              64
#> AKT3               1397             484             306             546
#> GAGE12F               0               0               0               0
#>         SAM80c6183220e6 SAM572f19794c96 SAM1f83ebd6be9b SAMe7e4f7c076a7
#> A1BG                 28              27               3              11
#> NAT2                  0               0               0               0
#> ADA                 934            1483              95             892
#> CDH2                171             170             422             549
#> AKT3                614             683             277            1298
#> GAGE12F               0               0               0               0
#>         SAMc6eff056c89a SAM5cfa1699bdb7 SAMda4d892fddc8 SAMe3d4266775a9
#> A1BG                 18             652              30              57
#> NAT2                  0               5               0               0
#> ADA                 652             215             738             339
#> CDH2                323            2089             358             126
#> AKT3               1012            2231             354             985
#> GAGE12F               0               0               0               0
```
