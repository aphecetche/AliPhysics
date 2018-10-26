# Quick assessment of AccxEff loss due to bus patch removal

The MCH Oncall sometimes has to decide whether or not to remove some bus patches from the configuration. The decision is a compromise between loss of tracking efficiency (by removing buspatches you might create "holes" in the acceptance) and operational burden (if the buspatches you are not removing are causing high busy times or worse crash the run(s)).

To make an informed decision one needs to be able to predict what loss of acceptance would result because of the removal of some buspatches.

That prediction can be made using some `AliMuonCompact*` classes.
You need four things to perform this prediction :

1. a file containing events from an ideal simulation (where everything is ok in the muon tracker)
2. a reference run number which serve as a basis to compute the AccxEff loss
3. a list of buspatches you're considering removing)
4. access to an OCDB

For instance, taking run 288806 as the reference, to compute the AccxEff loss that would result from the removal of buspatches 812 and 813 from the configuration, we'd do (assuming we get a local copy of 2018's OCDB) :

```
root[0] AliMuonCompactQuickAccEff qae("compacttreemaker.root")
root[1] AliMuonCompactManuStatus msref(288806,"local:///alice/data/2018/OCDB")
root[2] AliMuonCompactManuStatus ms = msref.RemoveBusPatches(std::vector<int>{812,813})
root[3] qae.ComputeLoss(msref,ms)
drop =    0.77 % +-  0.00 %
```

meaning the AccxEff would decrease by 0.77 % from its reference value in run 288806

```
root [4] msref.CountBad(ms.MANUOUTOFCONFIGMASK)
(int) 593
root[5] ms.CountBad(ms.MANUOUTOFCONFIGMASK)
(int) 615
```
