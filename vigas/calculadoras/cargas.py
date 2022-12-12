import numpy as np

class Carga:
    """Clase carga."""

    def __init__(self, tipo):
        """
        tipo = 0: Carga puntual
        tipo = 1: Carga distribuida
        tipo = 2: Momento cocentrado
        """
        self.tipo = tipo

    def tipo(self):
        if self.tipo == 0:
            print("Carga puntual")
        elif self.tipo == 1:
            print('Carga distribuida')
        elif self.tipo == 2:
            print('Momento concentrado')
        else:
            print('No definido')


class CargaPontual(Carga):
    '''Clase carga puntual'''
    def __init__(self, valor_carga=0, posicao=0):
        '''Carga puntual P.
        valor_carga: valor de la carga. Positivo hacia abajo.
        posicao: posicion de la carga respecto al extremo izquierdo del tramo.'''
        Carga.__init__(self, 0)
        self.valor_carga = valor_carga
        self.posicao = posicao
    
    def __str__(self):
        return 'Carga puntual\n   Valor= ' + str(self.valor_carga) + 'N' \
    + '\n   Posición, x= ' + str(self.posicao) + 'm'
    
    # Reacciones nodales equivalentes
    def Qf(self, longitude):
        '''Reacciones nodales equivalentes para una carga puntual.
        longitude: Longitud de la viga'''
        a = self.posicao
        b = longitude - a
        return self.valor_carga / longitude**2 * np.array([
                [b**2 / longitude * (3*a + b)],
                [a * b**2],
                [a**2 / longitude * (a + 3*b)],
                [-a**2 * b]
            ])
    
    # Fuerza cortante en una sección (viga sin apoyos)
    def FQ(self, posicao, longitude):
        '''Aporte a la fuerza cortante en una sección debido a una carga puntual,
        posicao: posición de la sección considerada respecto al extremo izquierdo
        longitude: longitud del tramo'''
        if self.posicao < posicao < longitude:
            return -self.valor_carga
        
        return 0
         
    #Momento flector en una sección (viga simplemente apoyada)
    def MF(self, posicao, longitude):
        '''Aporte al Momento flector en una sección debido a una carga puntual,
        posicao: posición de la sección considerada respecto al extremo izquierdo
        longitude: longitud del tramo'''
        if 0 <= posicao < self.posicao:
            return (1 - self.posicao/longitude) * self.valor_carga * posicao
        elif posicao <= longitude:
            return self.posicao * self.valor_carga * (1 - posicao/longitude)
        else:
            return 0



class CargaDistribuida(Carga):
    '''Clase carga distribuida'''
    def __init__(self, valor=0, posicao=0, longitude=0):
        '''Carga puntual P.
        valor: valor de la carga. Positivo hacia abajo.
        a: distancia entre el extremo izquierdo del tramo y el inicio de la carga.
        l: longitud de la carga distribuida'''
        Carga.__init__(self, 1)
        self.valor = valor
        self.inicio = posicao
        self.longitude = longitude
    
    def __str__(self):
        return 'Carga distribuida\n   Valor= ' + str(self.valor) + 'N/m'\
    ', ' + '\n   Inicio= ' + str(self.inicio) + 'm' + '\n   Longitud= ' + str(self.longitude) + 'm'
    
    def Qf(self, longitude):
        '''Reacciones nodales equivalentes para una carga
        unifomemente distribuida.
        longitude: longitud de la viga'''
        q = self.valor
        a = self.inicio
        b = longitude - self.inicio - self.longitude
        return q * longitude / 2 * np.array([
                [1 - a/longitude**4*(2*longitude**3 - 2*a**2*longitude + a**3) - b**3/longitude**4*(2*longitude - b)],
                [longitude/6*(1 - a**2/longitude**4*(6*longitude**2 - 8*a*longitude + 3*a**2) - b**3/longitude**4*(4*longitude - 3*b))],
                [1 - a**3/longitude**4*(2*longitude - a) - b/longitude**4*(2*longitude**3 - 2*b**2*longitude + a**3)],
                [-longitude/6*(1 - a**3/longitude**4*(4*longitude - 3*a) - b**2/longitude**4*(6*longitude**2 - 8*b*longitude + 3*b**2))]
            ])
            
    #Fuerza cortante en una sección (viga sin apoyos)
    def FQ(self, valor, longitude):
        '''Aporte a la fuerza cortante en una sección debido a la carga distribuida.
        valor: posición de la sección considerada respecto al extremo izquierdo
        longitude: Longitud del tramo'''
        if self.inicio <= valor < self.inicio + self.longitude:
            return -self.valor * (valor - self.inicio)
        elif valor <= longitude:
            return -self.valor * self.longitude
        else:
            return 0
    
    #Momento flector en una sección (viga simplemente apoyada)
    def MF(self, posicao, longitude):
        '''Aporte al momento flector en una sección debido a la carga distribuida.
        posicao: posición de la sección considerada respecto al extremo izquierdo
        longitude: Longitud del tramo'''
        V1 = self.valor*self.longitude/longitude*(longitude - self.inicio - self.longitude/2)
        V2 = self.valor*self.longitude - V1
        if 0 <= posicao < self.inicio:
            return V1 * posicao
        elif posicao <= self.inicio + self.longitude:
            return V1*posicao - 0.5*self.valor*(posicao-self.inicio)**2
        elif posicao <= longitude:
            return V2 * (longitude - posicao)
        else:
            return 0


class MomentoConcentrado(Carga):
    '''Clase momento concentrado'''
    def __init__(self, valor=0, posicao=0):
        '''Momento concentrado M.
        valor: valor del momento concentrado. Antihorario positivo
        posicao: posición del momento respecto al extremo izquierdo del tramo'''
        Carga.__init__(self, 2)
        self.valor = valor
        self.posicao = posicao
    
    def __str__(self):
        return 'Momento concentrado\n   Valor= ' + str(self.valor) + 'Nm' \
    + '\n   Posición, x= ' + str(self.posicao) + 'm'
    
    def Qf(self, longitude):
        '''Reacciones nodales equivalentes para un momento concetrado.
        longitude: longitud de la viga'''
        a = self.posicao
        b = longitude - a
        return self.valor / longitude**2 * np.array([
                [-6*a*b/longitude],
                [b*(b - 2*a)],
                [6*a*b/longitude],
                [a*(a - 2*b)]
            ])
    
    #Fuerza cortante en una sección (viga sin apoyos)
    def FQ(self, x, L):
        '''Aporte a la fuerza cortante en una sección debido a la carga distribuida.
        x: posición de la sección considerada respecto al extremo izquierdo'''
        return 0
    
    #Momento flector en una sección (viga simplemente apoyada)
    def MF(self, posicao, longitude):
        '''Aporte al momento flector en una sección debido a un momento concetrado,
        Estos valores corresponden al de una viga simplemente apoyada.
        posicao: posición de la sección considerada respecto al extremo izquierdo
        longitude: Longitud del tramo'''
        if 0 <= posicao < self.posicao:
            return self.valor / longitude * posicao
        elif self.posicao < posicao <= longitude:
            return self.valor * (posicao/longitude - 1)
        else:
            return 0