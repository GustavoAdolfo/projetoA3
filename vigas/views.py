import io
from functools import reduce
from pickle import LONG_BINPUT
import urllib
import base64
import numpy as np
import matplotlib.pyplot as plt
from urllib.error import HTTPError
from django.shortcuts import render
from django.http import HttpResponse, HttpResponseNotAllowed, HttpResponseBadRequest
from .calculadoras.bernoulli import VigaBernoulli
from .calculadoras.cargas import CargaPontual, MomentoConcentrado, CargaDistribuida

def calcular_matricial(request):
    if request.method != 'POST':
        return HttpResponseNotAllowed(request.method)

    nos = [int(no.replace('no_','')) for no in list(request.POST.keys()) if no.startswith('no_')]
    elasticidade = float(request.POST.get('inputElasticidade', 0))
    inercia = float(request.POST.get('inputInercia', 0))
    tamanho_viga = float(request.POST.get('inputTamanho', 0))
    apoio_esquerdo = int(request.POST.get('apoio_esquerdo', 0))
    apoio_direito = int(request.POST.get('apoio_direito', 0))

    pontos_na_viga = []
    cargas = []
    for no_ in nos:
        pontos_na_viga.append(VigaBernoulli(elasticidade, inercia, tamanho_viga))
        valor = float(request.POST.get(f'inputForca_{no_}', 0))
        posicao = float(request.POST.get(f'inputPosicao_{no_}', 0))
        dimensao = float(request.POST.get(f'inputDimensao_{no_}', 0))
        tipo_carga = int(request.POST.get(f'inputTipoCarga_{no_}', 0))
        if tipo_carga == 1:
            cargas.append([CargaPontual(valor, posicao)])
        elif tipo_carga == 2:
            cargas.append([CargaDistribuida(valor, posicao, dimensao)])
        elif tipo_carga == 3:
            cargas.append([MomentoConcentrado(valor, posicao)])
        else:
            pass

    qtde_pontos = len(pontos_na_viga)

    # Montando a matriz de rigidez global
    KG = np.zeros((2*(qtde_pontos+1), 2*(qtde_pontos+1)))
    for i in range(qtde_pontos):
        KG[2*i:2*i+4, 2*i:2*i+4] += pontos_na_viga[i].matriz_de_rigidez

    # Rigidez C
    C = np.amax(KG) * 1e4

    # Os graus de liberdade restringidos são:
    gdlRest = []

    # Em geral
    for i in range(qtde_pontos):
        gdlRest.append(2*i)
        
    # Extremo esquerdo
    if apoio_esquerdo == 0: # incorporação
        gdlRest.insert(1, 1)
    elif apoio_esquerdo == 1: # restrição de giro (permite rolagem vertical)
        gdlRest[0] = 1
    elif apoio_esquerdo == 3: # cantilever
        del gdlRest[0]
    else: # apoio de segundo grau
        pass

    # Extremo direito
    if apoio_direito == 0: # incorporação
        gdlRest.append(2*qtde_pontos)
        gdlRest.append(2*qtde_pontos + 1)
    elif apoio_direito == 1: # restrição de giro (permite rolagem vertical)
        gdlRest.append(2*qtde_pontos + 1)
    elif apoio_direito == 2: # apoio de segundo gênero (permite giro mas não deslocamento)
        gdlRest.append(2*qtde_pontos)
    else: # cantilever
        pass

    # Montagem da matriz de rigidez S pela abordagem de penalização
    S = KG
    for i in gdlRest:
        S[i,i] += C

    # Reacciones nodales equivalentes en cada tramo
    QF = [0]*qtde_pontos #para guardar los vectores de reacciones nodales equivalentes de cada tramo

    for i in range(qtde_pontos): # recorre todos los tramos
        for j in range(len(cargas[i])): # considera todas las cargas de cada tramo
            QF[i] += cargas[i][j].Qf(pontos_na_viga[i].longitude)

    # Montagem do vetor Qf para todos os gdl, incluídos os restringidos
    Qf = np.zeros((2*(qtde_pontos+1),1))
    for i in range(qtde_pontos):
        Qf[2*i:2*i+4,:] += QF[i]


    # Desplazamientos nodales
    d = -np.linalg.inv(S) @ Qf
    # Desplazamientos nodales por tramo
    u = []
    for i in range(qtde_pontos):
        u.append(d[2*i:2*i+4,:])

    #Fuerzas en cada tramo
    F = []
    for i in range(qtde_pontos):
        F.append(pontos_na_viga[i].matriz_de_rigidez @ u[i] + QF[i])
    
    # Reacciones
    r = d[gdlRest]
    R = -C * r

    #Número de secciones a tomar para los gráficos en cada tramo
    numS = 50
    Xt = [] #para guardar las x de cada tramo
    for i in range(qtde_pontos):
        Xt.append(np.linspace(0, pontos_na_viga[i].longitude, numS)) #Ubicación de las secciones

    Cortantes = []
    for i in range(qtde_pontos): #para cada tramo
        #Cortantes como vigas sin apoyo
        Q0 = np.zeros(numS)
        for j in range(len(cargas[i])): #considera todas las cargas de cada tramo
            m = 0 #para enumerar las secciones
            for x in Xt[i]: #recorre las secciones
                Q0[m] += cargas[i][j].FQ(x, pontos_na_viga[i].longitude)
                m += 1
        
        #Cortantes en el extremo, obtenidos del cálculo
        Q1 = F[i][0]
        
        #Momento total
        Cortantes.append(Q0 + Q1)
    
    #Máximos y mínimos valores de fuerza cortante (en cada tramo)
    maxCortante = [] #Cortantes máximos por cada tramo
    minCortante = [] #Cortantes mínimos por cada tramo
    XmaxQ= [] #ubicaciones de los máximos en cada tramo
    XminQ = [] #ubicaciones de los mínimos en cada tramo
    for i in range(qtde_pontos):
        maxQ = max(Cortantes[i]) #Máximo cortante
        minQ = min(Cortantes[i]) #Mínimo cortante
        maxCortante.append(maxQ)
        minCortante.append(minQ)
        indMaxQ = np.where(Cortantes[i] == maxQ )[0][0] #ubicación del máximo cortante
        indMinQ = np.where(Cortantes[i] == minQ )[0][0] #ubicación del mínimo cortante
        XmaxQ.append(Xt[i][indMaxQ])
        XminQ.append(Xt[i][indMinQ])

    Flectores = []
    for i in range(qtde_pontos): #para cada tramo            
        #Momentos como tramos simplemente apoyados
        M0 = np.zeros(numS)
        for j in range(len(cargas[i])): #considera todas las cargas de cada tramo
            m = 0 #para enumerar las secciones
            for x in Xt[i]: #recorre las secciones
                M0[m] += cargas[i][j].MF(x, pontos_na_viga[i].longitude)
                m += 1
        
        #Momentos debidos a los empotramientos o a la continuidad de la viga
        M1 = -F[i][1] + (F[i][3] + F[i][1]) / pontos_na_viga[i].longitude * Xt[i]
        
        #Momento total
        Flectores.append(M0 + M1)

    #Máximos y mínimos valores de momento flector (en cada tramo)
    maxFlector = [] #Flector máximo en cada tramo
    minFlector = [] #Flector mínimo en cada tramo
    XmaxF= [] #ubicaciones de los flectores máximos por tramo
    XminF = [] #ubicaciones de los mínimos flectores por tramo
    for i in range(qtde_pontos):
        maxF = max(Flectores[i]) #Máximo flector
        minF = min(Flectores[i]) #Mínimo flector
        maxFlector.append(maxF)
        minFlector.append(minF)
        indMaxF = np.where(Flectores[i] == maxF )[0][0] #ubicación del máximo flector
        indMinF = np.where(Flectores[i] == minF )[0][0] #ubicación del mínimo flector
        XmaxF.append(Xt[i][indMaxF])
        XminF.append(Xt[i][indMinF])

    ### 4.5 Diagrama de fuerza cortante

    #Valores de x para los gráficos
    X = []
    Lacum = 0
    for i in range(qtde_pontos):
        if i > 0:
            Lacum += pontos_na_viga[i-1].longitude
        Xprov = Xt[i] + Lacum
        Xlist = Xprov.tolist()
        X += Xlist

    #Valores de la fuerza cortante para los gráficos
    DFQ = []
    for i in range(qtde_pontos):
        #Valores para el DFQ tipo lista
        Corta = (Cortantes[i]/1000).tolist() #Pasamos a kN y convertimos en lista
        DFQ += Corta

    #Graf. principal de fuerza cortante
    plt.figure(1)
    plt.plot(X, DFQ)
    plt.title('Diagrama de Força Cortante', fontsize = 16)
    plt.xlabel('x [m]')
    plt.ylabel('Força cortante [kN]')
    plt.axhline(linewidth = 3)
    plt.xlim(0, tamanho_viga)
    plt.grid()

    colocarTextosQ(qtde_pontos, pontos_na_viga, XmaxQ, XminQ, tamanho_viga, maxCortante, minCortante)
    #Para sombrear el graf.
    Xgraf = [0] + X
    Xgraf.append(tamanho_viga)

    DFQgraf = [0] + DFQ
    DFQgraf.append(0)

    plt.fill(Xgraf, DFQgraf, 'b', alpha=0.3)

    #Divisores de tramos
    vertical = 0
    for i in range(qtde_pontos - 1):
        vertical += pontos_na_viga[i].longitude
        plt.axvline(vertical, color='black')

    # plt.show()
    fig = plt.gcf()
    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    buf.seek(0)
    string = base64.b64encode(buf.read())
    graf_forca_cortante = urllib.parse.quote(string)



    ### 4.4 Diagrama de momento flector
    #Valores del momento flector para los gráficos
    DMF = []
    for i in range(qtde_pontos):
        #Valores para el DMF tipo lista
        Flex = (Flectores[i]/1000).tolist() #Pasamos a kNm y convertimos en lista
        DMF += Flex

    #Graf. principal
    plt.figure(2)
    plt.plot(X, DMF)
    plt.title('Diagrama de Momento Fletor', fontsize = 16)
    plt.xlabel('x [m]')
    plt.ylabel('Momento fletor [kNm]')
    plt.gca().invert_yaxis() #invierte el eje y
    plt.axhline(linewidth = 3)
    plt.xlim(0, tamanho_viga)
    plt.grid()

    colocarTextosF(qtde_pontos, pontos_na_viga, XmaxF, XminF, tamanho_viga, maxFlector, minFlector)
            
    #Para sombrear el graf.
    Xgraf = [0] + X
    Xgraf.append(tamanho_viga)

    DMFgraf = [0] + DMF
    DMFgraf.append(0)

    plt.fill(Xgraf, DMFgraf, 'b', alpha=0.3)

    #Divisores de tramos
    vertical = 0
    for i in range(qtde_pontos - 1):
        vertical += pontos_na_viga[i].longitude
        plt.axvline(vertical, color='black')

    # plt.show()
    fig = plt.gcf()
    buf = io.BytesIO()
    fig.savefig(buf, format='png')
    buf.seek(0)
    string = base64.b64encode(buf.read())
    graf_momento_fletor = urllib.parse.quote(string)

    context = {
        'dados': {
            'forca_cortante': graf_forca_cortante,
            'momento_fletor': graf_momento_fletor
        }
    }

    return render(request, 'vigas/resultado_continua.html', context)



def continua(request):
    if request.method == 'GET':
        return render(request, 'vigas/continua.html')


#Función para colocar Textos de valores máximos y mínimos en flexión
def colocarTextosF(qtde_pontos, pontos_na_viga, XmaxF, XminF, tamanho_viga, maxFlector, minFlector):
    LacumM = 0
    for i in range(qtde_pontos):
        if i > 0:
            LacumM += pontos_na_viga[i-1].longitude
        ubicMax = LacumM + XmaxF[i]
        ubicMin = LacumM + XminF[i]
        if ubicMax == tamanho_viga:
            ubicMax = tamanho_viga - pontos_na_viga[i].longitude/2
        if ubicMin == tamanho_viga:
            ubicMin = tamanho_viga - pontos_na_viga[i].longitude/2
        plt.text(ubicMax*2.5, maxFlector[i]*0.00100, '$M_{max} = $' + \
                str(round(maxFlector[i]/1000,2)) + '$kNm, x= $' + str(round(XmaxF[i],2)) \
                + '$m$')
        plt.text(ubicMin*1.2, minFlector[i]*0.001, '$M_{min} = $' + \
                str(round(minFlector[i]/1000,2)) + '$kNm, x= $' + str(round(XminF[i],2)) \
                + '$m$')


#Textos para valores máximos y mínimos
def colocarTextosQ(qtde_pontos, pontos_na_viga, XmaxQ, XminQ, tamanho_viga, maxCortante, minCortante):
    LacumQ = 0
    for i in range(qtde_pontos):
        if i > 0:
            LacumQ += pontos_na_viga[i-1].longitude
        ubicMax = LacumQ + XmaxQ[i]
        ubicMin = LacumQ + XminQ[i]
        if ubicMax == tamanho_viga:
            ubicMax = tamanho_viga - pontos_na_viga[i].longitude/2
        if ubicMin == tamanho_viga:
            ubicMin = tamanho_viga - pontos_na_viga[i].longitude/2
        plt.text(ubicMax, maxCortante[i]*0.00108, '$Q_{max} = $' + \
                str(round(maxCortante[i]/1000,2)) + '$kN, x= $' + str(round(XmaxQ[i],2)) \
                + '$m$')
        plt.text(ubicMin, minCortante[i]*0.00108, '$Q_{min} = $' + \
                str(round(minCortante[i]/1000,2)) + '$kN, x= $' + str(round(XminQ[i],2)) \
                + '$m$')