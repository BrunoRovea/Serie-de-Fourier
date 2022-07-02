# Série de Fourier

Algoritmo que plota o resultado do somatório da série de Fourier, como dados de entrada são necessários:
* Período fundamental do sinal (T0)
* Número de termos da série (N)
* Número de harmõnicos que serão plotados no espectro em frequência (Ne)

O sinal é deescrito por duas funções tipo __handle__ (X1 e X2) dependentes da variável simbólica (t), sendo que a primeira função deescreve seu comportamento na primeira metade do período do sinal (de 0 a T0/2) e a segunda função tipo __handle__ do sinal deescreve seu comportamento na segunda metade do período (de T0/2 a T0).

Também é necessário especificar os intervalos de integração (I1 e I2) para obtenção dos coeficientes da série (an, bn, cn, dn).
