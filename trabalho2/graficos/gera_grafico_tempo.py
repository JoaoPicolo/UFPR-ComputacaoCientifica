import matplotlib.pyplot as plt

#Tempo
tempo_Otm_x1 = [10,32,50,64,100,128,200,256,300,400,512, 1000] #Otimizado
tempo_x = [10,32,50,64,100,128,200,256,300,400,512, 1000]
tempo_Otm_y1 = [0.000977,0.007275,0.020410,0.036865,0.122437,0.226074,0.835645,1.736898,2.800171,6.420736,12.693359, 139.777669]
tempo_y = [0.000488,0.010156,0.034473,0.068278,0.249023,0.515625,1.899072,3.897868,6.176514,14.390951,29.919678, 238.304688]

#L3
L3_x = [10,32,50,64,100,128,200,256,300,400,512, 1000] 
L3_Otm_x = [10,32,50,64,100,128,200,256,300,400,512, 1000] #Otimizado
L3_y = [22632.5375,2711.6476,1523.2869,1112.3925,1124.6049,614.2685,11714.8819,22261.8266,24724.9047,28046.9986,27263.1554, 24538.8219]
L3_Otm_y = [21117.4967,4244.3267,2553.7889,2036.8620,1524.2429,1548.7065,22279.8181,50832.1735,55246.5026,61112.8840,62719.7366, 40826.2735] #Otimizado

L2_x = [10,32,50,64,100,128,200,256,300,400,512, 1000] 
L2_Otm_x = [10,32,50,64,100,128,200,256,300,400,512, 1000] #Otimizado
L2_y = [0.2639,0.2600,0.2193,0.0756,0.0217,0.0217,0.1445,0.2155, 0.2387,0.2545,0.2365, 0.2537]
L2_Otm_y = [0.2141,0.2402,0.1985,0.0718,0.0204,0.0371,0.1441,0.2260,0.2342,0.2464,0.2332, 0.2474]

FLOPS_DP_x = [10,32,50,64,100,128,200,256,300,400,512, 1000] 
FLOPS_DP_Otm_x = [10,32,50,64,100,128,200,256,300,400,512, 1000] 
FLOPS_DP_y = [635.9844,1983.5897,2343.2688,2498.0728,2612.9111,2598.2295,2777.0825,2829.8673,2888.5652,2942.0343,2974.5811, 2723.1290]
FLOPS_DP_Otm_y = [734.3190,3204.9020,4152.6122,4890.7631,5366.8021,6288.0131,6537.4937,6352.9896,6403.7958,6680.5193,6961.4862, 4599.5458] 

FLOPS_AVX_x = [10,32,50,64,100,128,200,256,300,400,512, 1000] 
FLOPS_AVX_Otm_x = [10,32,50,64,100,128,200,256,300,400,512, 1000] 
FLOPS_AVX_y = [0,0,0,0,0,0,0,0,0,0,0, 0]
FLOPS_AVX_Otm_y = [647.6734,3168.5127,4024.4812, 4841.3758,5220.9263,6220.7103,6533.3787,6414.525,6270.5516,6657.8371,6925.9050, 4741.0297]

plt.title("FLOPS - Triangularização da Matriz U")
plt.ylabel('MFLOP/s', fontsize=10)
plt.xlabel('Tamanho de Entrada', fontsize=10)

# plt.yscale('log') 
plt.grid(True)
# plt.ylim(ymax=sorted(y)[-1]+100,ymin=sorted(y)[0]-1)

# Create scatter plot:
plt.plot(FLOPS_DP_x, FLOPS_DP_y, color="blue", marker="o",  linestyle="--" ,label="FLOPS_DP Não Otimizado")
plt.plot(FLOPS_DP_Otm_x, FLOPS_DP_Otm_y, color="red", marker="o",  linestyle="--" ,label="FLOPS_DP Otimizado")
plt.plot(FLOPS_AVX_x, FLOPS_AVX_y, color="navy", marker="o",  linestyle="--" ,label="FLOPS_AVX Não Otimizado")
plt.plot(FLOPS_AVX_Otm_x, FLOPS_AVX_Otm_y, color="darkred", marker="o",  linestyle="--" ,label="FLOPS_AVX Otimizado")

plt.legend(loc="lower right")
plt.xticks(L3_x,L3_x)
plt.show()