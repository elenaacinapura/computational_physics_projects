# Istruzioni per eseguire i programmi
Innanzitutto, eseguire da terminale i seguenti comandi, trovandosi all'interno della stessa directory di questo file README:
```
mkdir build
cd build
cmake ..
```
Poi seguire le istruzioni per i singoli esercizi, eseguendo i comandi riportati trovandosi sempre all'interno di `build`.

## Esercizio 1
```
make run-es1
```
## Esercizio 2
Per la prima parte, eseguire
```
make run-es2-1-numerov
```
per utilizzare il metodo di Numerov, oppure
```
make run-es2-1-diag
```
per utilizzare la diagonalizzazione.

Per la seconda parte, settare nel file `evolution.c` i parametri desiderati
- TYPE: 0 per simulazione semplice, 1 per simulatione pensata per graficare l'animazione dell'evoluzione temporale.
- POTENTIALTYPE: 0 per il primo potenziale (quello simmetrico), 1 per il secondo potenziale.
Poi, eseguire
```
make run-es2-2
```
per eseguire la simulazione. Per visualizzare i grafici:
- se si è eseguito con TYPE = 1, eseguire
```
make plot-es2-2-animation
```
- se si è eseguito con TYPE = 0, eseguire 
 ```
 make plot-es2-2-prob
 ```
 per visualizzare i grafici di p(t) e dello spettro.
## Esercizio 3
```
make run-es3
```

## Esercizio 4
```
make run-es4
```