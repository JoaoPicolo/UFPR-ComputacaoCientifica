import matplotlib.pyplot as plt

# X axis values:
x = [10,32,50,64,100,128,200,256,300,400,512, 1000] #Otimizado
x1 = [10,32,50,64,100,128,200,256,300,400,512, 1000]
# Y axis values:
y = [0.003906,0.053467,0.135742,0.220703,0.588989,0.949463,2.431641,3.825358,5.509949,9.820638,16.468343,64.594971]
y1 = [0.017090,0.201074,0.520410,0.822673,2.098022,3.406006,9.089697,14.244385,20.382080,36.448730, 61.771647, 243.407308]

plt.title("Tempo - Função Preenche Sistema Ajuste")
plt.ylabel('Tempo (ms)', fontsize=10)
plt.xlabel('Tamanho de Entrada', fontsize=10)

plt.yscale('log') 
plt.grid(True)
# plt.ylim(ymax=sorted(y)[-1]+100,ymin=sorted(y)[0]-1)

# Create scatter plot:
plt.plot(x, y, color="blue", marker="o",  linestyle="--" ,label="Otimizado")
# plt.plot(x1, y1, color="red", marker="o",  linestyle="--" ,label="Não Otimizado")

plt.legend(loc="upper left")
plt.xticks(x,x)
plt.show()