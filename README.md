# Charge Density Wave
Summary codes used in my CDW research \
Updating without a fixed schedule

## Phonon instability (q-dependent)
### Fermi Surface Nesting (FSN)
* Static Lindhard susceptibility 
```math
{\chi'}_\mathbf{q}=\sum_{\mathbf{k}}\frac{f(\varepsilon_\mathbf{k})-f(\varepsilon_{\mathbf{k}+\mathbf{q}})}{\varepsilon_{\mathbf{k}+\mathbf{q}}-\varepsilon_\mathbf{k}}
```
* Nesting function
```math
{\chi''}_\mathbf{k}=\sum_{\mathbf{q}}\delta\left(\varepsilon_\mathbf{k}\right)\delta\left(\varepsilon_{\mathbf{k}+\mathbf{q}}\right)
```
### Electron-Phonon Coupling (EPC)
* Moment-dependent electron-phonon coupling
```math
{\bar{g}}_\mathbf{q}=\sum_{\mathbf{k}}\left|g_{\mathbf{k},\mathbf{k}+\mathbf{q}}\right|^2
```
### FSN + EPC
* Generalized static electronic susceptibility
```math
\chi_\mathbf{q}=\sum_{\mathbf{k}}{\left|g_{\mathbf{k},\mathbf{k}+\mathbf{q}}\right|^2\frac{f(\varepsilon_\mathbf{k})-f(\varepsilon_{\mathbf{k}+\mathbf{q}})}{\varepsilon_{\mathbf{k}+\mathbf{q}}-\varepsilon_\mathbf{k}}},
```
## Scattered electrons (k-dependent)
Transfer sum over $\mathbf{q}$ to sum over $\mathbf{k}$ is from my institution, which should correspond to the scattered electrons (electron instability).
### Elastic scattered electrons
```math
{\chi'}_\mathbf{k}=\sum_{\mathbf{q}}\frac{f\left(\varepsilon_\mathbf{k}\right)-f\left(\varepsilon_{\mathbf{k}+\mathbf{q}}\right)}{\varepsilon_{\mathbf{k}+\mathbf{q}}-\varepsilon_\mathbf{k}}
```
```math
{\chi''}_\mathbf{k}=\sum_{\mathbf{q}}\delta\left(\varepsilon_\mathbf{k}\right)\delta\left(\varepsilon_{\mathbf{k}+\mathbf{q}}\right)
```
### Inelastic scattered electrons
```math
\chi_\mathbf{k}=\sum_{\mathbf{q}}{\left|g_{\mathbf{k},\mathbf{k}+\mathbf{q}}\right|^2\frac{f(\varepsilon_\mathbf{k})-f(\varepsilon_{\mathbf{k}+\mathbf{q}})}{\varepsilon_{\mathbf{k}+\mathbf{q}}-\varepsilon_\mathbf{k}}},
```


## Publications
If you use these CDW codes in your work, please consider citing:
>> Z. Wang, C. Chen, J. Mo, J. Zhou, K. P. Loh, and Y. P. Feng, _Decisive role of electron-phonon coupling for phonon and electron instabilities in transition metal dichalcogenides._
   [Phys. Rev. Research **5**, 013218 (2023)](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.5.013218)  
>> Z. Wang, J-Y. You, C. Chen, J. Mo, J. He, L. Zhang, J. Zhou, K. P. Loh, and Y. P. Feng, _Interplay of the charge density wave transition with topological and superconducting properties_
   [Nanoscale Horizons **8**, 1395 (2023)](https://pubs.rsc.org/en/content/articlelanding/2023/nh/d3nh00207a)
