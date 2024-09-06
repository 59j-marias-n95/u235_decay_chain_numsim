# Simulación numérica estocástica de la desintegración en cadena de U235
![Alt Text](https://github.com/59j-marias-n95/u235_decay_chain_numsim/blob/main/informe/figuras/cadenaluz.png)

Bienvenido al repositorio para el proyecto de investigación en el voluntariado subsecuente al curso de Introducción a la simulación numérica con GNU Octave propiciado por el Instituto de Investigación en Energía de la Universidad Nacional Autónoma de Honduras. 
En este proyecto se busca construir una simulación numérica del sistema de ecuaciones completo de la desintegración en cadena del radioisótopo Uranio-235.
Como elemento innovador en el estudio se plantea contemplar el sistema de ecuaciones como afectado por un ruido estocástico natural a la aleatoriedad del proceso. En vista de esta modificación se busca implementar métodos conocidos para resolver ecuaciones diferenciales acopladas de forma numérica con una variante estocástica. 
Los métodos a explorar son principalmente 2: el método Runge-Kutta y el método Euler-Maruyama, un tercer esquema de resolución es la simulación por Monte Carlo. 

## Plan de investigación
La resolución de este tipo de problemas se puede remontar prudentemente hasta Bateman quien estableció fórmulas generales para los coeficientes de integración de las ecuaciones diferenciales acopladas que modelan el proceso. Estas soluciones fueron abordadas por el método de transformada de Laplace que permite integrar las condiciones iniciales en la solución durante la misma ejecución del proceso y no en una etapa posterior. 
No obstante resolver un problema tan complejo computacionalmente como este mediante transformadas de Laplace resulta en un custo de cómputo elevado, no tanto por la búsqueda de la solución en el espacio de Laplace, pero por la transformación de las soluciones en dicho espacio hacia el espacio del tiempo otra vez. 
Otros métodos directos son consultados en la literatura y se encuentra que los esquemas matriciales son favorecidos. 

En este proyecto se planteó resolver el sistema de ecuaciones haciendo uso de los modelos de la teoría de Sistemas Dinámicos. Se planteó el sistema de ecuaciones como una ecuación matricial con una matriz de coeficientes, un vector dependiente del tiempo y su vector derivada. 
El sistema se resuelve buscando los eigenvalores de la matriz de coeficientes y computando sus eigenvectores para su inmediata linealización a través del vector de condiciones iniciales. 

Para establecer la viabilidad de este procedimiento se resuelve el sistema por transformadas de Laplace al no poder ejecutar las fórmulas de Bateman directamente a raíz de la complejida que subyace en la modificación de sus ecuaciones ante la consideración de tasas de ramificación. 
Laplace ejecutado de esta manera resulta ser agotar para la máquina en la ecuación número 15, por lo que la comparación de términos se contempla hasta la ecuación número 14 para indagar la solución número 15 en otro lenguaje o equipo de software. 
