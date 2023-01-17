######################################################################
##    Uso de HMM para añadir secuencias a una alineación múltiple   ##
######################################################################

# Borra todas las variables
rm(list=ls())
# Dice donde está el directorio de trabajo
getwd()
 


######################################################################
## Definir un conjunto de funciones que permitan:
## 1.- Obtener los datos de una alineación a partir del
##     contenido de un fichero de texto. Debe incluirse una
##     descripción del tipo de datos utilizado para contener toda la
##     información que se estime necesaria sobre la alineación.

##     El fichero que contenga la alineación debe tener el siguiente
##     formato ALM:

## -------------------------------------------------------------------
## #Nombre1
## [{L|-}]+
## #Nombre2
## [{L|-}]+
## ...
## -------------------------------------------------------------------
##     * La primera línea del fichero contiene el símbolo "#" seguido
##       de un nombre (sin espacios en blanco).
##     * La segunda línea contiene la primera de las secuencias que
##       componen la alineación: una sucesión de letras, sin espacios
##       entre ellas, que pueden contener el símbolo "-"
##     * Las restantes líneas, siempre en número par, tienen
##       sucesivamente el mismo formato que las dos primeras
##     * Todas las secuencias que componen la alineación son de la
##       misma longitud (el número de líneas y la longitud de las
##       secuencias puede variar de una alineación a otra). 




# ------------------------------------------------------------------------------
# Función que declara las secuencias dentro de un archivo de texto
obtenerAlineacion <- function(filepath) {
  # Función que lee el archivo
  con = file(filepath, "r")
  # Función que convierte el contenido del txt en una cadena en la que
  # cada elemento es una linea del txt.
  line = readLines(con)
  # Función que divide la cadena en dos. Una en la que cada elemento
  # es un nombre y otra en la que cada elemento es una secuencia.
  s <- split(line, 1:2)
  # Quitamos los corchetes a los nombres
  names <- substring(s[[1]], 2)
  sequences <- s[[2]]
  # Vector solución en la que almacenaremos las secuencias
  alignment <- c()
  # seq_along() es una función que crea un vector que contiene una
  # secuencia de números desde 1 hasta la longitud del argumento.
  for (i in seq_along(sequences)){
    # coge la secuencia i y la separa en una lista de caracteres
    sequence = strsplit(sequences[i], "")[1]
    alignment <- c(alignment, sequence)
  }
  # Se asigna su nombre a cada secuencia
  names(alignment) <- names
  return(alignment)
}
# ------------------------------------------------------------------------------
alineacion1 <- obtenerAlineacion('alineacionMultiple18.txt')
alineacion1
# ---------


# ------------------------------------------------------------------------------
# Función que dada una alineación devuelve la secuencia de regiones de inserción
# y coincidencia
obtenerRegiones <- function(alineacion){
  # Vector de solución
  regiones <- c()
  # Número de la m
  contador_m <- 1
  # Número de la i
  contador_i <- 0
  # Por cada una de las posiciones de una secuencia
  for (i in seq_along(alineacion[[1]])){
    # Se inicia un contador para letras y guiones
    letter_counter <- 0
    hyphen_counter <- 0
    # Por cada una de las secuencias
    for (j in seq_along(alineacion)){
      # Si en la posición i de la secuencia j hay un guion se aumenta el
      # contador de guiones
      if(alineacion[[j]][i]=="-"){
        hyphen_counter <- hyphen_counter + 1
      # Si hay una letra se aumenta el contador de letras
      }else{
        letter_counter <- letter_counter + 1
      }
    }
    # Si hay más letras se añade una M
    if (letter_counter > hyphen_counter){
      regiones[i] <- paste("M", contador_m, sep="")
      contador_m <- contador_m + 1
    # En caso contrario se añade una I
    }else{
      # Si la región anterior ya es una I, la posición actual es la misma que la
      # anterior
      if (substr(regiones[i-1], 1, 1) == "I"){
        regiones[i] <- regiones[i-1]
      # En caso contrario es una I nueva
      }else{
        regiones[i] <- paste("I", contador_i, sep="")
        contador_i <- contador_i + 1
      }
    }
  }
  return(regiones)
}
# ------------------------------------------------------------------------------
regiones1 <- obtenerRegiones(alineacion1)
regiones1
# ---------



# ------------------------------------------------------------------------------
# Función que nos devuelve todas las letras que aparecen en la alineación
obtenerAlfabeto <- function(alineacion){
  alfabeto <- c()
  for (i in alineacion){
    for (j in i){
      # Si la letra no es un guión y no está en el alfabeto, se añade al alfabeto
      if (j != "-" & j %in% alfabeto == FALSE){
        alfabeto <- c(alfabeto, j)
      }
    }
  }
  # Devuelve el alfabeto ordenado
  return(sort(alfabeto))
}

# ------------------------------------------------------------------------------
alfabeto1 <- obtenerAlfabeto(alineacion1)
alfabeto1
# ---------





## 2.- Construir un modelo oculto de Markov a partir de los
##     datos de una alineación, a partir de la información obtenida
##     sobre una alineación y por tanto, según el tipo de datos
##     descrito en el punto anterior. Debe incluirse una descripción
##     del tipo de datos utilizado para contener toda la información
##     que se estime necesaria sobre el modelo.




# ------------------------------------------------------------------------------
# Función que dada una alineación devuelve el vector de probabilidades de empezar
# por cada región
vectorInicial <- function(alineacion){
  # Obtenemos las regiones correspondientes a la alineación
  regiones <- obtenerRegiones(alineacion)
  # Vector en el que recopilar las veces en la que cada región empieza por una 
  # letra
  vector_contador <- c()
  for (i in seq_along(alineacion[[1]])){
    vector_contador[i] <- 0
  }
  # Si la secuencia i presenta su primera letra en la posición j entonces se
  # suma 1 en la posición j del vector contador y se pasa a la siguiente secuencia
  for (i in seq_along(alineacion)){
    for (j in seq_along(alineacion[[1]])){
      if (alineacion[[i]][j] != "-"){
        vector_contador[j] <- 1 + vector_contador[j]
        break
      }
    }
  }
  # Vector de probabilidades, coge cada elemento del vector_contador y se dividen
  # entre la suma total de sus elementos.
  vector_inicial <-c()
  for (i in seq_along(regiones)){
    vector_inicial[i] <- vector_contador[i]/sum(vector_contador)
  }
  return(vector_inicial)
}
# ------------------------------------------------------------------------------
vector_inicial1 <- vectorInicial(alineacion1)
vector_inicial1
# ---------




# ------------------------------------------------------------------------------
# Función que dada una alineación devuelve la matriz de transición en la que
# aparecen las probabilidades de transición de una región a otra.

# Se ejemplifica lo que hace la función con este ejemplo en el que se calcula la
# primera fila de la matriz de probabilidades:
# Vector de posiciones: 3 2 2 3
# Alineación
# "E" "-" "R" "E" "-" "N" "-" "-"
# "-" "-" "-" "-" "-" "K" "-" "-"
# "E" "D" "I" "V" "W" "M" "R" "T"
# "E" "I" "L" "-" "-" "M" "-" "-"
# "E" "-" "L" "K" "-" "M" "-" "-"
# Regiones
# "M1" "I0" "M2" "M3" "I1" "M4" "I2" "I2"
# En la primera secuencia, la primera letra es E y la siguiente es R en la
# posición 3 que está en la siguiente región: M2
# En la segunda secuencia, hay una primera letra pero no una segunda
# En la tercera secuencia, la primera letra es E y la siguiente es D en la
# posición 2 que está en la siguiente región: I0
# En la cuarta secuencia, la primera letra es E y la siguiente es I en la
# posición 2 que está en la siguiente región: I0
# En la quinta secuencia, la primera letra es E y la siguiente es L en la
# posición 3 que está en la siguiente región: M2

# Las probabilidades son:
# Hay 2/5 de que la siguiente sea I0
# Hay 2/5 de que la siguiente sea M2

# En la matriz:    
#    M1  I0  M2   M3        I1        M4 I2 I2
# M1  0 0.5 0.5 0.00 0.0000000 0.0000000  0  0

matrizTransicion <- function(alineacion){
  # Longitud de una secuencia, coincide con longitud de la región
  longitud <- length(alineacion[[1]])
  # Lista en la que almacenaremos las posiciones de las letras
  lista <- list()
  for (i in 1:longitud){
    lista[[i]] <-0
  }
  # Bucle desde 1 hasta la longitud de la secuencia menos 1. No necesitamos 
  # averiguar la siguiente letra a la situada en la última posición
  for (i in 1:(longitud-1)){
    # Vector auxiliar en el que almacenar las posiciones
    vector <- c();
    # Contador auxiliar para añadir elementos al vector
    counter <- 1;
    # Bucle que va de 1 hasta la longitud de la alineación
    for (j in seq_along(alineacion)){
      aux <- i + 1
      # Sí en la secuencia hay un guión, se pasa a la siguiente
      if (alineacion[[j]][i] == "-"){
        next
        # Sí hay una letra:
      }else{
        # Se busca la posición de la siguiente letra. Por eso el bucle va de 
        # i+1 hasta la longitud de la secuencia
        for (k in aux:longitud){
          # Si lo que hay en la posición k de la alineación j no es un guión
          if (alineacion[[j]][k] != "-"){
            # Se añade al vector la posición correspondiente a la siguiente letra
            vector[counter] <- k
            counter <- counter+1
            break
          }
        }
      }
    }
    # Se almacena el vector en la lista de posiciones
    lista[[i]] <- vector
  }
  # Obtenemos las regiones correspondientes a la alineación
  regiones <- obtenerRegiones(alineacion)
  # Creamos una matriz vacía con tantas filas y columnas como regiones
  matriz_transicion <- matrix(0, ncol = length(regiones), nrow = length(regiones)) 
  # Le asignamos a las filas y a las columnas los nombres de las regiones
  colnames(matriz_transicion) <- regiones 
  rownames(matriz_transicion) <- regiones
  # Rellenamos la matriz con las probabilidades
  for (i in 1:(longitud)){
    total <- length(lista[[i]])
    for (j in 1:(longitud)){
      matriz_transicion[i,lista[[i]][j]] <- matriz_transicion[i,lista[[i]][j]] + 1/total
    }
  }
  return(matriz_transicion)
}
# ------------------------------------------------------------------------------
matriz_transicion1 <- matrizTransicion(alineacion1)
matriz_transicion1 
# ---------


# ------------------------------------------------------------------------------
# Función que dada una alineación devuelve las probabilidades con la que cada uno
# de las regiones presentan unos observables u otros.
matrizDeObservacion <- function(alineacion){
  # Obtenemos las regiones correspondientes a la alineación
  regiones <- obtenerRegiones(alineacion)
  # Obtenemos el alfabeto
  alfabeto <- obtenerAlfabeto(alineacion)
  # Creamos una matriz vacía con tantas filas y columnas como regiones
  matriz_de_observacion <- matrix(0, ncol = length(alfabeto), nrow = length(regiones)) 
  # Le asignamos a las filas y a las columnas los nombres de las regiones
  colnames(matriz_de_observacion) <- alfabeto 
  rownames(matriz_de_observacion) <- regiones
  # Bucle que va desde i hasta el número de regiones
  for (i in seq_along(regiones)){
    # Vector de ceros del tamaño del alfabeto
    vector_aux <- replicate(length(alfabeto), 0)
    # Bucle que va desde i hasta el número de alineaciones
    for (j in seq_along(alineacion)){
      # Bucle que va desde k hasta el número de letras en el alfabeto
      for (k in seq_along(alfabeto)){
        # Si la letra en la posición i de la secuencia j es la misma que la letra
        # del alfabeto en la posición k entonces se suma uno en la posición k del
        # vector auxiliar
        if (alineacion[[j]][i] == alfabeto[k]){
          vector_aux[k] <- vector_aux[k] + 1
        }
      }
    }
    # Se sustituye cada fila de la matriz final por el vector auxiliar dividido
    # entre la suma de sus componentes
    for (j in seq_along(alfabeto)){
      matriz_de_observacion[i, j] <-  vector_aux[j]/sum(vector_aux)
    }
  }
  return(matriz_de_observacion)
}
# ------------------------------------------------------------------------------
matriz_de_observacion1 <- matrizDeObservacion(alineacion1)
matriz_de_observacion1
# ---------


## 3.- Añadir una secuencia a la alineación, utilizando para
##     ello la sucesión más probable de regiones dada la secuencia,
##     según el modelo anterior. Hay que tener en cuenta que no hay
##     límite para la longitud de las secuencias de la alineación y de
##     la secuencia a añadir.

##     Nota: es suficiente con que la función proporcione dicha
##     sucesión de regiones.


viterbi <- function (alineacion, secuencia){
  # Modelo de Markov
  vector_inicial <- vectorInicial(alineacion)
  matriz_transicion <- matrizTransicion(alineacion)
  matriz_de_observacion <- matrizDeObservacion(alineacion)
  # Regiones de la alinecion
  regiones <- obtenerRegiones(alineacion)
  # Longitud de la secuencia
  longitud_secuencia <- length(secuencia)
  # Matriz con número de filas longitud de la secuencia menos 1, número de
  # columnas longitud de las regiones. El comando dimnames le da nombre a filas
  # y columnas. Primero a las filas y luego a las columnas
  pr <- matrix(nrow = longitud_secuencia-1, ncol = length(regiones),
               dimnames = list(2:longitud_secuencia, regiones))
  # nu_1(e_j)
  nu <- vector_inicial*matriz_de_observacion[,secuencia[1]]
  for (k in 2:longitud_secuencia) {
    # mt * nu_k-1
    aux <- matriz_transicion * nu
    # apply(X, MARGIN, FUN) aplica una función a todos los elementos de una matriz
    # X: Una matriz
    # MARGIN: Aplicar en: 1 son filas y 2 son colummnas.
    # FUN: La función que aplicaremos a la matriz X en filas o columnas
    # nu_{k}(e_j) = mo_{e_js_k}*max(mt_{e,e_j}v_{k-1}(e))
    nu <- matriz_de_observacion[,secuencia[k]] * apply(aux, 2, max)
    # pr_k(e_j) = argmax(mt_{e,e_j}v_{k-1}(e))
    pr[as.character(k),] <- apply(aux, 2, which.max)
  }
  # q_t = argmax(v_t(e_j))
  q <- c()
  q[longitud_secuencia] <- regiones[which.max(nu)]
  # para k desde t-1 hasta 1
  for (k in (longitud_secuencia-1):1) {
    # q_k = pr_k+1(q_k+1)
    q[k] <- regiones[pr[as.character(k+1),q[k+1]]]
  }
  return(q)
}

secuencia_referencia1 <- c("E", "E", "N")
viterbi(alineacion1, secuencia_referencia1)




viterbiLogaritmo <- function (alineacion, secuencia) {
  # Modelo de Markov y aplicarle el logaritmo
  vector_inicial <- vectorInicial(alineacion)
  vector_inicial_log <- log(vector_inicial)
  matriz_transicion <- matrizTransicion(alineacion)
  matriz_transicion_log <- log(matriz_transicion)
  matriz_de_observacion <- matrizDeObservacion(alineacion)
  matriz_de_observacion_log <- log(matriz_de_observacion)
  # Longitud de la secuencia
  longitud_secuencia <- length(secuencia)
  # Regiones de la alinecion
  regiones <- obtenerRegiones(alineacion)
  # Matriz con número de filas longitud de la secuencia menos 1, número de
  # columnas longitud de las regiones. El comando dimnames le da nombre a filas
  # y columnas. Primero a las filas y luego a las columnas
  pr <- matrix(nrow = longitud_secuencia-1, ncol = length(regiones),
               dimnames = list(2:longitud_secuencia, regiones))
  # nu_1(e_j)
  nu <- matriz_de_observacion_log[,secuencia[1]] + vector_inicial_log
  #
  for (k in 2:longitud_secuencia){
    # v_{k-1}(e) + log(mt_{e,e_j}) 
    aux <- nu + matriz_transicion_log
    # apply(X, MARGIN, FUN) aplica una función a todos los elementos de una matriz
    # X: Una matriz
    # MARGIN: Aplicar en: 1 son filas y 2 son colummnas.
    # FUN: La función que aplicaremos a la matriz X en filas o columnas
    # nu_{k}(e_j) = log(mo_{e_js_k}) + max(log(mt_{e,e_j}) + v_{k-1}(e))
    nu <- matriz_de_observacion_log[,secuencia[k]] + apply(aux, 2, max)
    # pr_k(e_j) = argmax(log(mt_{e,e_j}) + v_{k-1}(e))
    pr[as.character(k),] <- apply(aux, 2, which.max)
  }
  # q_t = argmax(v_t(e_j))
  q <- c()
  # La función which devuelve la posición si el argumento de la función es TRUE
  q[longitud_secuencia] <- regiones[which.max(nu)]
  # para k desde t-1 hasta 1
  for (k in (longitud_secuencia-1):1){
    # q_k = pr_k+1(q_k+1)
    q[k] <- regiones[pr[as.character(k+1),q[k+1]]]
  }
  return(q)
}

secuencia_referencia1 <- c("E", "E", "N")
qq <- viterbiLogaritmo(alineacion1, secuencia_referencia1)
qq



######################################################################
## (1) Descargar de la web de Uniprot, en un único fichero con formato
##     FASTA, las secuencias del conjunto de proteínas asignado. A
##     partir de este, obtener (se recomienda utilizar MegaX, pero
##     puede usarse cualquier otra herramienta incluyendo la
##     referencia correspondiente) una alineación múltiple de las
##     secuencias del conjunto. Guardar la alineación en un fichero de
##     texto y, si es necesario, editarlo para adaptarlo al formato
##     descrito en el apartado anterior.
## (2) Obtener un modelo oculto de Markov asociado a la alineación
##     anterior utilizando la correspondiente función definida en el
##     primer apartado


# ---------
alineacion2 <- obtenerAlineacion('alineacion_multiple_proteinas.txt')
alineacion2
# ---------
regiones2 <- obtenerRegiones(alineacion2)
regiones2
# ---------
alfabeto2 <- obtenerAlfabeto(alineacion2)
alfabeto2
# ---------
vector_inicial2 <- vectorInicial(alineacion2)
vector_inicial2
# ---------
matriz_transicion2 <- matrizTransicion(alineacion2)
matriz_transicion2 
# ---------
matriz_de_observacion2 <- matrizDeObservacion(alineacion2)
matriz_de_observacion2
# ---------

## (3) Descargar de la web de Uniprot, en un fichero con formato FASTA
##     independiente del obtenido en el punto anterior, la secuencia
##     de la proteína de referencia asignada.

library (seqinr)

datos <- getSequence(read.fasta("data.fasta"))
sec_1 <- toupper(datos[[1]])
sec_2 <- toupper(datos[[2]])
sec_3 <- toupper(datos[[3]])
sec_4 <- toupper(datos[[4]])
sec_5 <- toupper(datos[[5]])
sec_6 <- toupper(datos[[6]])
secuencia_referencia2 <- toupper(datos[[7]])
secuencia_referencia2

## (4) Añadir la proteína de referencia a la alineación del conjunto
##     de proteínas utilizando la correspondiente función definida en
##     el apartado anterior.
######################################################################

qq2 <- viterbi(alineacion2, secuencia_referencia2)
qq2
qq2_log <- viterbiLogaritmo(alineacion2, secuencia_referencia2)
qq2_log

qq2_sec_1 <- viterbiLogaritmo(alineacion2, sec_1)
qq2_sec_1

qq2_sec_2 <- viterbiLogaritmo(alineacion2, sec_2)
qq2_sec_2

qq2_sec_3 <- viterbiLogaritmo(alineacion2, sec_3)
qq2_sec_3

qq2_sec_4 <- viterbiLogaritmo(alineacion2, sec_4)
qq2_sec_4

qq2_sec_5 <- viterbiLogaritmo(alineacion2, sec_5)
qq2_sec_5

qq2_sec_6 <- viterbiLogaritmo(alineacion2, sec_6)
qq2_sec_6




