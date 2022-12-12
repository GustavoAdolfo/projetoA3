import numpy as np
import matplotlib.pyplot as plt


class VigaBernoulli:
    '''"Viga de Bernoulli" - definida pelas suas características de \
        módulo de elasticidade, inercia e longitude; com esses dados \
            se obtém a matriz de rigidez do elemento.

    Definimos un tramo de viga.
    E: Módulo de elasticidad
    I: Inercia de la sección transversal
    L: Longitud del tramo'''
    def __init__(self, elasticidade, inercia, longitude):
        """Definição de uma seção de viga.

        ATRIBUTOS:
            self.E: Módulo de elasticidad
            self.I: Inercia de la sección transversal
            self.L: Longitud del tramo
            self.k: matriz de rigidez del tramo
        """
        self.elasticidade = elasticidade
        self.inercia = inercia
        self.longitude = longitude

        # Matriz de rigidez do elemento
        self.matriz_de_rigidez = elasticidade * inercia / longitude**3 * np.array([
                [12., 6*longitude, -12, 6*longitude],
                [6*longitude, 4*longitude**2, -6*longitude, 2*longitude**2],
                [-12, -6*longitude, 12, -6*longitude],
                [6*longitude, 2*longitude**2, -6*longitude, 4*longitude**2]
            ])