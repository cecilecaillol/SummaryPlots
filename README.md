Plot of limits as a function of ma. You can choose the model and if the limit is on BR(h->aa) times BR(a->mumu)^2 or times BR(a->tautau)^2 (--br tau or --br mu). The plots are the same in type 1 and type 2 and do not depend on tan beta. In type-3 and type-4, mmbb analysis is drawn for tan beta = 1, 2 and 4

```shell
python Plot_allModels_allMa.py --model 4 --br tau
python Plot_allModels_allMa.py --model 3 --br tau
python Plot_allModels_allMa.py --model 1 --br mu
python Plot_allModels_allMa.py --model 1 --br tau
```

Plot limit for a given mass as a function of tan beta (interesting in type 3 and type 4) for mmbb and mmtt (so mass between 20 and 62.5 GeV)

```shell
python Plot_mmtt_mmbb_allTanBeta.py --model 3 --tanbeta 1 --ma 40
python Plot_mmtt_mmbb_allTanBeta.py --model 4 --tanbeta 1 --ma 40
```

Plot all branching fractions as a function of the mass, depending on the model and tan beta. The mass is used only in the print-out but does not change the plots.

```shell
python Plot_BR.py --model 1 --tanbeta 5 --ma 40 #no dependence on tan beta in type 1
python Plot_BR.py --model 2 --tanbeta 5 --ma 40
python Plot_BR.py --model 3 --tanbeta 5 --ma 40
python Plot_BR.py --model 4 --tanbeta 5 --ma 40
```

Plot BR(aa-->XXXX) as a function of tan beta for a given mass in all types of 2HDM+S

```shell
python Plot_BR_allTanBeta.py --channel mmmm --ma 2
python Plot_BR_allTanBeta.py --channel tttt --ma 7
python Plot_BR_allTanBeta.py --channel tttt --ma 15
python Plot_BR_allTanBeta.py --channel mmbb --ma 40
python Plot_BR_allTanBeta.py --channel mmtt --ma 40
```

Plot limit on BR(h->aa). Inputs are the type of model and tan beta (tan beta not used in type-1).

```shell
python Plot_BRaa.py --model 1
python Plot_BRaa.py --model 2 --tanbeta 2
python Plot_BRaa.py --model 3 --tanbeta 5
python Plot_BRaa.py --model 4 --tanbeta 0.5
```

Plot bbA analyses comparison

```shell
python Plot_bbA.py
```

