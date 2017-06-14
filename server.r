library(shiny)
require(mvtnorm)
require(igraph)
require(visNetwork)
require(plotrix)

shinyServer(
  function(input, output){
    seed <- reactive(input$seed)
    dataselect <- reactive(input$dataselect)
    alphaselect <- reactive(input$alphaselect)
    deltax <- reactive(input$x)
    deltay <- reactive(input$y)
    letop <- reactive(input$top)
    
    # I4CS
    
     
######### Les fonctions : rareness(), detect_top() ##########
    
    figure <- reactive({
      source <- seed()
      set.seed(source)
      # set.seed(91)
      
      data_n <- NULL
      n <- 0
      while(n<47){
        datacube <- c(runif(n = 1, min=1, max = 9), 
                      runif(n = 1, min=1, max = 9))
        rayon <- (datacube[1]-5)*(datacube[1]-5) + (datacube[2]-5)*(datacube[2]-5)
        datacreux <- as.numeric( (4<rayon) && (rayon<16) )
        if(datacreux != 0){
          n <- n+1
          data_n <- rbind(data_n, datacube)
        }
      }
      data_n <- rbind(data_n, c(5,5))
      data_n <- rbind(data_n, c(1, 9))
      data_n <- rbind(data_n, c(4, 8))
      
      # numero de la donnee : l'individu
      num <- dataselect()
      dx <- deltax()
      dy <- deltay()
      data_n[num, 1] <- data_n[num, 1]+dx
      data_n[num, 2] <- data_n[num, 2]+dy
      
      data <- data_n
      data
    })
    
    
    obs <- function(data){
      n <- dim(data)[1]   # nb de donnees
      p <- dim(data)[2]   # nb d attributs ou dimension
      # Calcul des p matrices de distances par attribut 
      # p matrices n * n
      # calculee sur les lignes de data (par c(2) dans le apply())
      lesDistances <- array(apply(data, c(2), 
                                  function(x){ as.matrix( 
                                    dist(x, diag = TRUE, upper = TRUE, method= "euclidean"))}),
                            c(n, n, p))
      # Calcul des observations
      # Mettre les elements diagonaux a -1 pour eviter les egalites
      # pour chacune des p matrices de distances (par c(3) dans le apply())
      lesDistDiag <- array(apply(lesDistances, c(3), 
                                 function(x){ y<-as.matrix(x); diag(y) <- -1; return(y) }), 
                           c(n, n, p))
      # calcul des rangs : les observations
      # Attention les rangs obtenus seront en colonnes (et non en lignes)
      lesObservations <- array(apply(lesDistDiag, c(2,3), 
                                     function(x){rank(x, ties.method="first")}), 
                               c(n, n, p))
      return(lesObservations)
    }
      
    rareness <- function(data, lesObservations, alpha, threshold){
      n <- dim(data)[1]   # nb de donnees
      p <- dim(data)[2]   # nb d attributs ou dimension
        
      # isolement pour l'individu i dans le voisinage j
      # l'isolement = les p rangs de l'individu i vu depuis l'observateur k
      isolement <- matrix( rep(0, n*p), c(n, p) )
      # Pour chaque donnees i
      for(i in 1:n){
        # n voisinages
        # Calcul des voisinages (indicateur en colonne pour chaque attribut)
        lesVoisinages <- array(apply(lesObservations[,i,], c(1,2), 
                                     function(x){((1<x)&&(x<=alpha))}), 
                               c(n,p))
        # Pour chaque voisinage (attribut j)
        for(j in 1:p){
          # p attributs
          # les rangs individu i pour son voisinage j (en colonne)
          # avec alpha-1 element en colonne pour le voisinage j
          # et p colonnes (p attributs)
          lesRangs <- lesObservations[i, which(lesVoisinages[,j]==TRUE),]
          # tronquer les rangs par threshold
          lesRangs <- pmin(lesRangs, threshold)
          
          #lesRangs[lesRangs==threshold] <- 0
          
          lesRangs <- matrix(lesRangs, alpha-1, p)
          # les coefficients isolement de i dans son voisinage j
          # moyenne des voisins pour chacun des p attributs
          lesCoeffIsolement <- apply(lesRangs, 2, mean)
          # Prendre l'isolement pour l'attribut isolant le plus i
          # (on reste dans le voisinage j)
          isolement[i,j] <- max(lesCoeffIsolement/(alpha))
        }
      }
      return(isolement)
    }
    
    total_rare <- function(data, lesObservations, alpha){
      n <- dim(data)[1]   # nb de donnees
      p <- dim(data)[2]   # nb d attributs ou dimension
      isolAttribut <- matrix(rareness(data, lesObservations, alpha, alpha), n, p)
      coeffIsol <- apply(isolAttribut, 1, max)
      tot <- sum(coeffIsol == 1)
      #if(sum(coeffIsol<1) ==  0){
      #  signif <- 1
      #}else{
      #  signif <- sum(coeffIsol[coeffIsol<1])/sum(coeffIsol<1)
      #}
      #total <- cbind(tot, signif)
      return(tot)
    }
      
    detect_alpha <- function(data, lesObservations, top, debut, fin){
      top <- top-1
      n <- dim(data)[1]   # nb de donnees
      p <- dim(data)[2]   # nb d attributs ou dimension
      print(paste("Recherche alpha optimal entre", debut, "et", fin))

      # heuristique pour trouver alpha
      totdebut <- total_rare(data, lesObservations, debut)
      totfin <- total_rare(data, lesObservations, fin)
      milieu <- trunc((debut+fin)/2)
      totmilieu <- total_rare(data, lesObservations, milieu)
      #print(paste("debut", debut, "(", totdebut[1], totdebut[2], "), ",
      #           "milieu", milieu, "(", totmilieu[1], totmilieu[2], "), ",
      #           "fin", fin, "(", totfin[1], totfin[2], ")"
      #           ))
      print(paste("", debut, "(", totdebut, ")", " < ", 
                      milieu, "(", totmilieu, ")", " < ", 
                      fin, "(", totfin, ")"))
      while(milieu > debut){
        #if( (totmilieu[1] > top) && (totmilieu[2] > 0.95) ){
        if( (totmilieu > top) ){
            debut <- milieu
          totdebut <- totmilieu
        }else{
          fin <- milieu
          totfin <- totmilieu
        }
        milieu <- trunc((debut+fin)/2)
        totmilieu <- total_rare(data, lesObservations, milieu)
        print(paste("", debut, "(", totdebut, ")", " < ", 
                    milieu, "(", totmilieu, ")", " < ", 
                    fin, "(", totfin, ")"))
      }
      alpha <- debut
      if((totdebut >= top)&&(totfin != 0)){
        alpha <- fin
      }
      print(paste("  Alpha optimal pour le top", top+1, ":", alpha))
      return(alpha)
    }
    
    affiche_rareness <- function(data, lesObservations, alpha, seuil){
      n <- dim(data)[1]   # nb de donnees
      p <- dim(data)[2]   # nb d attributs ou dimension
      isolAttribut <- matrix(rareness(data, lesObservations, alpha, seuil), n, p)
      return(isolAttribut)
    }
    
    
    
    
############### Figure : Data Set ##########
    
    output$dataset <- renderPlot({
      data <- figure()
      # print("Afficher les donnees : ")
      # print(data)

      n <- dim(data)[1]   # nb de donnees
      p <- dim(data)[2]   # nb d attributs ou dimension
      
      fig <- data[,1:2]
      colnames(fig) <- c("Attribute 1", "Attribute 2")
      par(mfrow=c(1,2))
      
      plot(fig, 
           xlim=c(0, 10), 
           ylim=c(0, 10),
           cex=1.3, pch=3, lwd=2)
      title(main=paste("2D Dataset ( n =", n, ")"))
      labels <- c("C","A","B")
      j <- 1
      for(i in 48:50){
        text(data[i,1], data[i,2], label=labels[j], cex= 1.5, pos=4)
        j <- j+1
        points(fig[i,1], fig[i,2], lty="solid", pch=3, cex=1.5, lwd=4, col="red")
      }
      draw.ellipse(5, 5, a=2, b=2, lty="dashed")
      draw.ellipse(5, 5, a=4, b=4, lty="dashed")
      
      plot(fig, 
           xlim=c(0, 10), 
           ylim=c(0, 10),
           xaxt='n', yaxt='n',
           cex=1.3, pch=3, lwd=2)
      title(main=paste("Obs: Observer of the Dataset"))
      # numero de la donnee : l'individu
      num <- dataselect()
      x <- data[num, 1]
      y <- data[num, 2]
      points(x, y, lty="solid", pch=3, cex=1.5, lwd=3, col="blue")
      text(x, y, label=" Obs", cex= 1.5, pos=4, col="blue")
      segments(x, -2, x, y, lty="dashed", lwd=4, col = "blue")
      lab_X <- c(0, round(x, 2), 10)
      axis(side=1, at=c(0, x, 10), cex.axis = 1.2,
           labels=lab_X, 
           lwd=1, lwd.ticks=2, col.ticks="blue")
      segments(-2, y, x, y, lty="dashed", lwd=4, col = "blue")
      lab_Y <- c(0, round(y, 2), 10)
      axis(side=2, at=c(0, y, 10), cex.axis = 1.2,
           labels=lab_Y, 
           lwd=1, lwd.ticks=2, col.ticks="blue")
      
      })

    
    
    
    
############### Figure ranking : Observations ##########
    
    output$ranking <- renderPlot({
      data <- figure()
      n <- dim(data)[1]   # nb de donnees
      p <- dim(data)[2]   # nb d attributs ou dimension
      
      # les rangs
      lesObservations <- obs(data)
      # numero de la donnee : l'individu
      num <- dataselect()
      fig <- data[,1:2]
      colnames(fig) <- c("Attribute 1", "Attribute 2")

      par(mfrow=c(1,2))
      
      plot(fig, 
           xlim=c(0, 10), 
           ylim=c(0, 10),
           cex=1.3, pch=3, lwd=2,
           xlab="Ranks of A, B and C (Ranking from Obs)",
           ylab="Attribute 2",
           xaxt='n')
      axis(1, labels=FALSE, col = "red")
      title(main=paste("Observation with Attribute 1"))
      labels <- c("C","A","B")
      j <- 1
      for(i in 48:50){
        text(data[i,1], data[i,2], label=labels[j], cex= 1.5, pos=4)
        j <- j+1
        points(fig[i,1], fig[i,2], lty="solid", pch=3, cex=1.3, lwd=3, col="red")
      }
      x <- data[num,1]
      y <- data[num,2]
      points(x, y, lty="solid", pch=3, cex=1.5, lwd=4, col="blue")
      segments(x, -2, x, y, lty="dashed", lwd=4, col = "red")
      text(data[num,1], data[num,2], label=" Obs", cex= 1.5, pos=4, col="blue")
      for(i in 48:50){
        x <- data[i,1]
        y <- data[i,2]
        segments(x, -2, x, y, lty="dashed", lwd=2, col = "red")
        points(x, y, lty="solid", pch=3, cex=1.3, lwd=2, col="red")
      }
      lab_C <- lesObservations[48, num, 1]
      lab_A <- lesObservations[49, num, 1]
      lab_B <- lesObservations[50, num, 1]
      lab_X <- lesObservations[num, num, 1]
      axis(side=1, at=data[48,1],
           cex.axis = 1.2,
           labels=lab_C,
           lwd=1, lwd.ticks=2, col.ticks="red")
      axis(side=1, at=data[49,1],
           cex.axis = 1.2,
           labels=lab_A,
           lwd=1, lwd.ticks=2, col.ticks="red")
      axis(side=1, at=data[50,1],
           cex.axis = 1.2,
           labels=lab_B,
           lwd=1, lwd.ticks=2, col.ticks="red")
      axis(side=1, at=c(data[num,1]), 
           cex.axis = 1.2,
           labels=lab_X, 
           lwd=1, lwd.ticks=4, col.ticks="red")
      
      
      plot(fig, 
           xlim=c(0, 10), 
           ylim=c(0, 10),
           cex=1.3, pch=3, lwd=2,
           ylab="Ranks of A, B and C (Ranking from Obs)",
           xlab="Attribute 1",
           yaxt='n')
      axis(2, labels=FALSE, col = "red")
      title(main=paste("Observation with Attribute 2"))
      labels <- c("C","A","B")
      j <- 1
      for(i in 48:50){
        text(data[i,1], data[i,2], label=labels[j], cex= 1.5, pos=4)
        j <- j+1
        points(fig[i,1], fig[i,2], lty="solid", pch=3, cex=1.3, lwd=3, col="red")
      }
      x <- data[num,1]
      y <- data[num,2]
      points(x, y, lty="solid", pch=3, cex=1.5, lwd=4, col="blue")
      segments(-2, y, x, y, lty="dashed", lwd=4, col = "red")
      text(data[num,1], data[num,2], label=" Obs", cex= 1.5, pos=4, col="blue")
      for(i in 48:50){
        x <- data[i,1]
        y <- data[i,2]
        segments(-2, y, x, y, lty="dashed", lwd=2, col = "red")
        points(x, y, lty="solid", pch=3, cex=1.3, lwd=2, col="red")
      }
      lab_C = lesObservations[48, num, 2]
      lab_A = lesObservations[49, num, 2]
      lab_B = lesObservations[50, num, 2]
      lab_X = lesObservations[num, num, 2]
      axis(side=2, at=data[48,2], 
           cex.axis = 1.2,
           labels=lab_C, 
           lwd=1,lwd.ticks=2,col.ticks="red")
      axis(side=2, at=data[49,2], 
           cex.axis = 1.2,
           labels=lab_A, 
           lwd=1,lwd.ticks=2,col.ticks="red")
      axis(side=2, at=data[50,2], 
           cex.axis = 1.2,
           labels=lab_B, 
           lwd=1,lwd.ticks=2,col.ticks="red")
      axis(side=2, at=c(data[num,2]), 
           cex.axis = 1.2,
           labels=lab_X, 
           lwd=1,lwd.ticks=4,col.ticks="red")
      
    })
    
    
    
    
    
############### Figure : Neighborhoods ##########
    
    output$neighbors <- renderPlot({
      data <- figure()
      n <- dim(data)[1]   # nb de donnees
      p <- dim(data)[2]   # nb d attributs ou dimension
      
      # les rangs
      lesObservations <- obs(data)
      # les alpha plus proches voisins
      alpha <- alphaselect()
      # numero de la donnee : l'individu
      num <- dataselect()
      # voisinage pour l'individu i
      # Calcul des voisinages (indicateur en colonne pour chaque attribut)
      lesVoisinages <- array(apply(lesObservations[,num,], c(1,2), 
                                   function(x){((1<x)&&(x<=alpha))}), 
                             c(n,p))
      fig <- data[,1:2]

      par(mfrow=c(1,2))
      
      colnames(fig) <- c("Neighborhood with Attribute 1", "Attribute 2")
      plot(fig, 
           xlim=c(0, 10), 
           ylim=c(0, 10),
           cex=1.3, pch=3, lwd=2)
      title(main=paste("k Nearest Neighbors with Attribute 1 ( k =",alpha,")"))
      labels <- c("C","A","B")
      j <- 1
      for(i in 48:50){
        text(data[i,1], data[i,2], label=labels[j], cex= 1.5, pos=4)
        j <- j+1
        points(fig[i,1], fig[i,2], lty="solid", pch=3, cex=1.3, lwd=3, col="red")
      }
      x <- data[num,1]
      y <- data[num,2]
      points(x, y, lty="solid", pch=3, cex=1.5, lwd=4, col="blue")
      # segments(x, -2, x, y, lty="dashed", lwd=4, col = "blue")
      text(data[num,1], data[num,2], label=" Obs", cex= 1.5, pos=4, col="blue")
      
      neighbors_x <-data[which(lesVoisinages[,1]==TRUE), ]
      for(i in 1:(alpha-1)){
        x <- neighbors_x[i,1]
        y <- neighbors_x[i,2]
        # segments(x, -2, x, y, lty="dashed", lwd=2, col = "blue")
        points(x, y, lty="solid", pch=3, cex=1.5, lwd=3, col="blue")
      }
      x <- max(neighbors_x[,1])
      segments(x, -5, x, 10.5, lwd=2, col = "blue")
      x <- min(neighbors_x[,1])
      segments(x, -5, x, 10.5, lwd=2, col = "blue")
      
      
      colnames(fig) <- c("Attribute 1", "Neighborhood with Attribute 2")
      plot(fig, 
           xlim=c(0, 10), 
           ylim=c(0, 10),
           cex=1.3, pch=3, lwd=2)
      title(main=paste("k Nearest Neighbors with Attribute 2 ( k =",alpha,")"))
      labels <- c("C","A","B")
      j <- 1
      for(i in 48:50){
        text(data[i,1], data[i,2], label=labels[j], cex= 1.5, pos=4)
        j <- j+1
        points(fig[i,1], fig[i,2], lty="solid", pch=3, cex=1.3, lwd=3, col="red")
      }
      x <- data[num,1]
      y <- data[num,2]
      points(x, y, lty="solid", pch=3, cex=1.5, lwd=4, col="blue")
      # segments(-2, y, x, y, lty="dashed", lwd=4, col = "blue")
      text(data[num,1], data[num,2], label=" Obs", cex= 1.5, pos=4, col="blue")
      
      neighbors_y <-data[which(lesVoisinages[,2]==TRUE), ]
      for(i in 1:(alpha-1)){
        x <- neighbors_y[i,1]
        y <- neighbors_y[i,2]
        # segments(-2, y, x, y, lty="dashed", lwd=2, col = "blue")
        points(x, y, lty="solid", pch=3, cex=1.5, lwd=3, col="blue")
      }
      y <- max(neighbors_y[,2])
      segments(-5, y, 10.5, y, lwd=2, col = "blue")
      y <- min(neighbors_y[,2])
      segments(-5, y, 10.5, y, lwd=2, col = "blue")
      
    })
    

    
    
############### Figure : Neighbors are observers ##########
    
    output$neighbor <- renderPlot({
      data <- figure()
      n <- dim(data)[1]   # nb de donnees
      p <- dim(data)[2]   # nb d attributs ou dimension
      
      # les rangs
      lesObservations <- obs(data)
      # numero de la donnee : l'individu
      num <- dataselect()
      # les alpha plus proches voisins
      alpha <- alphaselect()
      fig <- data[,1:2]
      colnames(fig) <- c("x", "y")
      lesVoisinages <- array(apply(lesObservations[,num,], c(1,2), 
                                   function(x){((1<x)&&(x<=alpha))}), 
                             c(n,p))
      
      par(mfrow=c(1,2))
      
      plot(fig, 
           xlim=c(0, 10), 
           ylim=c(0, 10),
           cex=1.3, pch=3, lwd=2,
           xlab="Rank of Obs from N",
           ylab="Rank of Obs from N",
           xaxt='n', yaxt='n')
      title(main=paste("Rank of Obs by N with Attribute 1"))
      labels <- c("C","A","B")
      j <- 1
      for(i in 48:50){
        text(data[i,1], data[i,2], label=labels[j], cex= 1.5, pos=4)
        j <- j+1
        points(fig[i,1], fig[i,2], lty="solid", pch=3, cex=1.3, lwd=3, col="red")
      }
      
      neighbors_x <-data[which(lesVoisinages[,1]==TRUE), ]
      for(i in 1:(alpha-1)){
        x <- neighbors_x[i,1]
        y <- neighbors_x[i,2]
        # segments(x, -2, x, y, lty="dashed", lwd=2, col = "blue")
        points(x, y, lty="solid", pch=3, cex=1.5, lwd=3, col="blue")
      }
      x <- max(neighbors_x[,1])
      segments(x, -5, x, 10.5, lty="dashed", lwd=1, col = "blue")
      x <- min(neighbors_x[,1])
      segments(x, -5, x, 10.5, lty="dashed", lwd=1, col = "blue")
      
      extremity <- c(which(lesObservations[,num,1] == alpha))
      
      # Pour l'observateur
      x <- data[num,1]
      y <- data[num,2]
      segments(x, -2, x, y, lty="dashed", lwd=4, col = "red")
      points(x, y, lty="solid", pch=3, cex=1.5, lwd=3, col="blue")
      text(data[num,1], data[num,2], label=" Obs", cex= 1.5, pos=4, col="blue")
      lab_X <- lesObservations[num, extremity, 1]
      if(lab_X > alpha){
        lab_X <- paste("", lab_X, "( threshold to", alpha, ")")
      }
      axis(side=1, at=c(data[num,1]), cex.axis = 1.2,
           labels=lab_X, 
           lwd=1, lwd.ticks=4, col.ticks="red")

      # Pour le voisin
      neighbor_x <-data[extremity, ]
      x <- neighbor_x[1]
      y <- neighbor_x[2]
      segments(x, -2, x, y, lty="dashed", lwd=4, col = "red")
      points(x, y, lty="solid", pch=3, cex=1.5, lwd=4, col="red")
      text(x, y, label="N", cex= 1.5, pos=4, col="black")
      # lab_X <- lesObservations[extremity, extremity, 1]
      lab_X <- "1"
      axis(side=1, at=c(x), cex.axis = 1.2,
           labels=lab_X, 
           lwd=1, lwd.ticks=4, col.ticks="red")

      
           
      plot(fig, 
           xlim=c(0, 10), 
           ylim=c(0, 10),
           cex=1.3, pch=3, lwd=2,
           ylab="Rank of Obs from N",
           xlab="Rank of Obs from N",
           xaxt='n', yaxt='n')
      title(main=paste("Rank of Obs by N with attribute 2"))
      labels <- c("C","A","B")
      j <- 1
      for(i in 48:50){
        text(data[i,1], data[i,2], label=labels[j], cex= 1.5, pos=4)
        j <- j+1
        points(fig[i,1], fig[i,2], lty="solid", pch=3, cex=1.3, lwd=3, col="red")
      }
      
      neighbors_x <-data[which(lesVoisinages[,1]==TRUE), ]
      for(i in 1:(alpha-1)){
        x <- neighbors_x[i,1]
        y <- neighbors_x[i,2]
        # segments(x, -2, x, y, lty="dashed", lwd=2, col = "blue")
        points(x, y, lty="solid", pch=3, cex=1.5, lwd=3, col="blue")
      }
      x <- max(neighbors_x[,1])
      segments(x, -5, x, 10.5, lty="dashed", lwd=1, col = "blue")
      x <- min(neighbors_x[,1])
      segments(x, -5, x, 10.5, lty="dashed", lwd=1, col = "blue")
      
      extremity <- c(which(lesObservations[,num,1] == alpha))
      
      # Pour l'observateur
      x <- data[num,1]
      y <- data[num,2]
      segments(-2, y, x, y, lty="dashed", lwd=4, col = "red")
      points(x, y, lty="solid", pch=3, cex=1.5, lwd=3, col="blue")
      text(data[num,1], data[num,2], label=" Obs", cex= 1.5, pos=4, col="blue")
      lab_X <- lesObservations[num, extremity, 2]
      if(lab_X > alpha){
        lab_X <- paste("", lab_X, "( threshold to", alpha, ")")
      }
      axis(side=2, at=c(data[num,2]), cex.axis = 1.2,
           labels=lab_X, 
           lwd=1, lwd.ticks=4, col.ticks="red")

      # Pour le voisin
      neighbor_x <-data[extremity, ]
      x <- neighbor_x[1]
      y <- neighbor_x[2]
      segments(-2, y, x, y, lty="dashed", lwd=4, col = "red")
      points(x, y, lty="solid", pch=3, cex=1.5, lwd=4, col="red")
      text(x, y, label="N", cex= 1.5, pos=4, col="black")
      lab_X <- "1"
      axis(side=2, at=c(y), cex.axis = 1.2,
           labels=lab_X, 
           lwd=1, lwd.ticks=4, col.ticks="red")
      
    })
    
    
     

    
############### Figure : Rareness ##########
    
    output$rareness <- renderPlot({
      data <- figure()
      n <- dim(data)[1]   # nb de donnees
      p <- dim(data)[2]   # nb d attributs ou dimension
      
      # les rangs
      lesObservations <- obs(data)
      # numero de la donnee : l'individu
      num <- dataselect()
      # les alpha plus proches voisins
      alpha <- alphaselect()
      # voisinage pour l'individu i
      # Calcul des voisinages (indicateur en colonne pour chaque attribut)
      lesVoisinages <- array(apply(lesObservations[,num,], c(1,2), 
                                   function(x){((1<x)&&(x<=alpha))}), 
                             c(n,p))
      fig <- data[,1:2]
     
      par(mfrow=c(1,2))
      
      plot(fig, 
           xlim=c(0, 10), 
           ylim=c(0, 10),
           xlab="Normalized Mean of thresholded Ranks of Obs",
           ylab="Normalized Mean of thresholded Ranks of Obs",
           xaxt='n', yaxt='n',
           cex=1.3, pch=3, lwd=2)
      title(main=paste("Neighborhood with Attribute 1 ( k =",alpha,")"))
      labels <- c("C","A","B")
      j <- 1
      for(i in 48:50){
        text(data[i,1], data[i,2], label=labels[j], cex= 1.5, pos=4)
        j <- j+1
        points(fig[i,1], fig[i,2], lty="solid", pch=3, cex=1.3, lwd=3, col="red")
      }
      
      x <- data[num,1]
      y <- data[num,2]
      points(x, y, lty="solid", pch=3, cex=1.5, lwd=4, col="blue")
      segments(x, -2, x, y, lty="dashed", lwd=4, col = "red")
      segments(-2, y, x, y, lty="dashed", lwd=4, col = "red")
      text(data[num,1], data[num,2], label=" Obs", cex= 1.5, pos=4, col="blue")
      
      lesRangs <- lesObservations[num, which(lesVoisinages[,1]==TRUE),]
      # tronquer les rangs par threshold=alpha
      lesRangs <- pmin(lesRangs, alpha)
      # matrice des rangs
      lesRangs <- matrix(lesRangs, alpha-1, p)
      # moyenne des voisins pour chacun des p attributs
      lesCoeffIsolement <- apply(lesRangs, 2, mean)
      lab <- round(lesCoeffIsolement[1]/alpha, 2)
      axis(side=1, at=data[num,1], cex.axis = 1.2,
           labels=lab,
           lwd=1, lwd.ticks=4, col.ticks="red")
      lab <- round(lesCoeffIsolement[2]/alpha, 2)
      axis(side=2, at=data[num,2], cex.axis = 1.2,
           labels=lab,
           lwd=1, lwd.ticks=4, col.ticks="red")
      
      neighbors_x <-data[which(lesVoisinages[,1]==TRUE), ]
      x <- max(neighbors_x[,1])
      segments(x, -5, x, 10.5, lty="dashed", lwd=3, col = "blue")
      x <- min(neighbors_x[,1])
      segments(x, -5, x, 10.5, lty="dashed", lwd=3, col = "blue")
      for(i in 1:(alpha-1)){
        x <- neighbors_x[i,1]
        y <- neighbors_x[i,2]
        segments(x, -2, x, y, lty="dashed", lwd=1, col = "red")
        segments(-2, y, x, y, lty="dashed", lwd=1, col = "red")
        points(x, y, lty="solid", pch=3, cex=1.3, lwd=2, col="blue")
      }
      lab <- rep("", (alpha-1))
      axis(side=1, at=neighbors_x[,1], cex.axis = 1.2,
           labels=lab,
           lwd=1, lwd.ticks=1, col.ticks="red")
      axis(side=2, at=neighbors_x[,2], cex.axis = 1.2,
           labels=lab,
           lwd=1, lwd.ticks=1, col.ticks="red")

      
      
      
      plot(fig, 
           xlim=c(0, 10), 
           ylim=c(0, 10),
           xlab="Normalized Mean of thresholded Ranks of Obs",
           ylab="Normalized Mean of Thresholded Ranks of Obs",
           xaxt='n', yaxt='n',
           cex=1.3, pch=3, lwd=2)
      title(main=paste("Neighborhood with Attribute 2 ( k =",alpha,")"))
      labels <- c("C","A","B")
      j <- 1
      for(i in 48:50){
        text(data[i,1], data[i,2], label=labels[j], cex= 1.5, pos=4)
        j <- j+1
        points(fig[i,1], fig[i,2], lty="solid", pch=3, cex=1.3, lwd=3, col="red")
      }
      
      x <- data[num,1]
      y <- data[num,2]
      points(x, y, lty="solid", pch=3, cex=1.5, lwd=4, col="blue")
      segments(x, -2, x, y, lty="dashed", lwd=4, col = "red")
      segments(-2, y, x, y, lty="dashed", lwd=4, col = "red")
      text(data[num,1], data[num,2], label=" Obs", cex= 1.5, pos=4, col="blue")
      
      lesRangs <- lesObservations[num, which(lesVoisinages[,2]==TRUE),]
      # tronquer les rangs par threshold=alpha
      lesRangs <- pmin(lesRangs, alpha)
      # matrice des rangs
      lesRangs <- matrix(lesRangs, alpha-1, p)
      # moyenne des voisins pour chacun des p attributs
      lesCoeffIsolement <- apply(lesRangs, 2, mean)
      lab <- round(lesCoeffIsolement[1]/alpha, 2)
      axis(side=1, at=data[num,1], cex.axis = 1.2,
           labels=lab,
           lwd=1, lwd.ticks=4, col.ticks="red")
      lab <- round(lesCoeffIsolement[2]/alpha, 2)
      axis(side=2, at=data[num,2], cex.axis = 1.2,
           labels=lab,
           lwd=1, lwd.ticks=4, col.ticks="red")
      
      neighbors_y <-data[which(lesVoisinages[,2]==TRUE), ]
      y <- max(neighbors_y[,2])
      segments(-5, y, 10.5, y, lty="dashed", lwd=3, col = "blue")
      y <- min(neighbors_y[,2])
      segments(-5, y, 10.5, y, lty="dashed", lwd=3, col = "blue")
      for(i in 1:(alpha-1)){
        x <- neighbors_y[i,1]
        y <- neighbors_y[i,2]
        segments(x, -2, x, y, lty="dashed", lwd=1, col = "red")
        segments(-2, y, x, y, lty="dashed", lwd=1, col = "red")
        points(x, y, lty="solid", pch=3, cex=1.3, lwd=2, col="blue")
      }
      lab <- rep("", (alpha-1))
      axis(side=1, at=neighbors_y[,1], cex.axis = 1.2,
           labels=lab,
           lwd=1, lwd.ticks=1, col.ticks="red")
      axis(side=2, at=neighbors_y[,2], cex.axis = 1.2,
           labels=lab,
           lwd=1, lwd.ticks=1, col.ticks="red")

    })
    

    
    
    
     

############### Figure : Special Cases ##########
    
    output$cases <- renderPlot({
      data <- figure()
      # print("Afficher les donnees : ")
      # print(data)
      
      n <- dim(data)[1]   # nb de donnees
      p <- dim(data)[2]   # nb d attributs ou dimension
      # numero de la donnee : l'individu
      num <- dataselect()
      # les alpha plus proches voisins
      alpha <- alphaselect()

      lesObservations <- obs(data)
      isolAttribut <- affiche_rareness(data, lesObservations, alpha,  alpha)
      coeffIsol <- apply(isolAttribut, 1, max)
      specialCases <- which(coeffIsol==1)
      letop <- sum(coeffIsol[specialCases])
      
      par(mfrow=c(1,2))
      
      fig <- data[,1:2]
      colnames(fig) <- c("Attribute 1", "Attribute 2")
      plot(fig, 
           xlim=c(0, 10), 
           ylim=c(0, 10),
           xlab="Two Means of thresholded Ranks of Obs",
           ylab="Two Means of Thresholded Ranks of Obs",
           xaxt='n', yaxt='n',
           cex=1.3, pch=3, lwd=2)
      title(main=paste("rareness(Obs) =", round(coeffIsol[num],2) ))
      labels <- c("C","A","B")
      j <- 1
      for(i in 48:50){
        text(data[i,1], data[i,2], label=labels[j], cex= 1.5, pos=4)
        j <- j+1
        points(fig[i,1], fig[i,2], lty="solid", pch=3, cex=1.3, lwd=3, col="red")
      }
      x <- data[num,1]
      y <- data[num,2]
      if(coeffIsol[num] == 1){
        points(x, y, lty="solid", pch=1, cex=2, lwd=3, col="red")
      }else{
        points(x, y, lty="solid", pch=3, cex=1.5, lwd=5, col="blue")
      }
      segments(x, -2, x, y, lty="solid", lwd=1, col = "blue")
      segments(-2, y, x, y, lty="solid", lwd=1, col = "blue")
      text(data[num,1], data[num,2], label=" Obs", cex= 1.5, pos=4, col="blue")
      lesVoisinages <- array(apply(lesObservations[,num,], c(1,2), 
                                   function(x){((1<x)&&(x<=alpha))}), 
                             c(n,p))
      neighbors_x <-data[which(lesVoisinages[,1]==TRUE), ]
      for(i in 1:(alpha-1)){
        x <- neighbors_x[i,1]
        y <- neighbors_x[i,2]
        points(x, y, lty="solid", pch=3, cex=1.3, lwd=3, col="blue")
        # segments(x, -2, x, y, lty="solid", lwd=0.5, col = "blue")
      }
      x <- max(neighbors_x[,1])
      segments(x, -5, x, 10.5, lty="dashed", lwd=3, col = "blue")
      x <- min(neighbors_x[,1])
      segments(x, -5, x, 10.5, lty="dashed", lwd=3, col = "blue")
      
      neighbors_y <-data[which(lesVoisinages[,2]==TRUE), ]
      for(i in 1:(alpha-1)){
        x <- neighbors_y[i,1]
        y <- neighbors_y[i,2]
        points(x, y, lty="solid", pch=3, cex=1.3, lwd=3, col="blue")
        # segments(-2, y, x, y, lty="solid", lwd=0.5, col = "blue")
      }
      y <- max(neighbors_y[,2])
      segments(-5, y, 10.5, y, lty="dashed", lwd=3, col = "blue")
      y <- min(neighbors_y[,2])
      segments(-5, y, 10.5, y, lty="dashed", lwd=3, col = "blue")
      
      
      lesRangs <- lesObservations[num, which(lesVoisinages[,1]==TRUE),]
      # tronquer les rangs par threshold=alpha
      lesRangs <- pmin(lesRangs, alpha)
      # matrice des rangs
      lesRangs <- matrix(lesRangs, alpha-1, p)
      # moyenne des voisins pour chacun des p attributs
      lesCoeffIsolement <- apply(lesRangs, 2, mean)
      lab_Nx_x <- round(lesCoeffIsolement[1]/alpha, 2)
      lab_Nx_y <- round(lesCoeffIsolement[2]/alpha, 2)
      
      lesRangs <- lesObservations[num, which(lesVoisinages[,2]==TRUE),]
      # tronquer les rangs par threshold=alpha
      lesRangs <- pmin(lesRangs, alpha)
      # matrice des rangs
      lesRangs <- matrix(lesRangs, alpha-1, p)
      # moyenne des voisins pour chacun des p attributs
      lesCoeffIsolement <- apply(lesRangs, 2, mean)
      lab_Ny_x <- round(lesCoeffIsolement[1]/alpha, 2)
      lab_Ny_y <- round(lesCoeffIsolement[2]/alpha, 2)
      
      axis(side=1, at=data[num,1], cex.axis = 1.2,
           labels=paste("(", lab_Nx_x, ",", lab_Ny_x, ")"),
           lwd=1, lwd.ticks=4, col.ticks="blue")
      axis(side=2, at=data[num,2], cex.axis = 1.2,
           labels=paste("(", lab_Nx_y, ",", lab_Ny_y, ")"),,
           lwd=1, lwd.ticks=4, col.ticks="blue")
      
                
      plot(fig, 
           xlim=c(0, 10), 
           ylim=c(0, 10),
           cex=1.3, pch=3, lwd=2)
      title(main=paste("",  letop, "Special Cases with rareness = 1 ( k = ", alpha, ")"))
      # title(sub=paste("Three special cases obtained using 5 nearest neighbours"))
      labels <- c("C","A","B")
      j <- 1
      for(i in 48:50){
        text(data[i,1], data[i,2], label=labels[j], cex= 1.5, pos=4)
        j <- j+1
        points(fig[i,1], fig[i,2], lty="solid", pch=3, cex=1.3, lwd=2, col="red")
      }
      for(i in specialCases){
        points(data[i,1], data[i,2], lty="solid", col="red", lwd=3, pch=1, cex=2)
        # text(data[i,1], data[i,2], label="1", cex= 1.5, pos=1)
      }
    })
    
    
    
    
    
############### Figure alpha : Size of Neighborhoods ##########
    
    output$alpha <- renderPlot({
      data <- figure()
      n <- dim(data)[1]   # nb de donnees
      p <- dim(data)[2]   # nb d attributs ou dimension
      lesObservations <- obs(data)
      
      alphaVal <- 2:n
      lesCas <- rep(0, length(alphaVal))
      i <- 1
      for(alpha in alphaVal){
        isolAttribut <- affiche_rareness(data, lesObservations, alpha, alpha)
        coeffIsol <- apply(isolAttribut, 1, max)
        lesCas[i] <- sum(coeffIsol[coeffIsol == 1])
        print(paste("Nb cas = ", lesCas[i], "pour alpha = ", alpha))
        i <- i+1
      }
      sauvegarde <- cbind(alphaVal, lesCas)
      # write.table(x=sauvegarde, file="courbe_data", 
      #            row.names=FALSE, col.names=FALSE, sep = "\t", quote=FALSE)
      # sauvegarde <- read.table(file="courbe_data", header = FALSE, sep = "\t")
      
      plot(sauvegarde, type="l", 
           xlab="Size of neighbourhoods ( k )", 
           ylab="Number of Special Cases",
           main="Number of Special Cases and Size of Neighbourhoods",
           lwd=3)
      
      x <- alphaselect()
      y <- lesCas[x-1]
      segments(x, -1, x, y, lty="dashed", lwd=2, col = "red")
      axis(1, at = c(x), labels = paste("", x), col.axis="red", 
           tick=TRUE, col="red", lwd=2)
      segments(-1, y, x, y, lty="dashed", lwd=2, col = "red")
      axis(2, at = c(y), labels = paste("", y), col.axis="red", 
           tick=TRUE, col="red", lwd=2)
      
      })
    
    
    
      
      
######### L'onglet Read Me (A COMPLETER) ##########
      
      output$readme <- renderUI({
        str1 <- h3("Detection of Special Cases")
        str2 <- "Concept of observer in a set of multidimensional data."
        str3 <- "Simulation with 50 data in 2D space
                 (each data has two attributes x and y)."
        str4 <- "The attributes are processed separately."
        str5 <- "The special cases are data that is inconsistent with its neighborhoods."
        HTML(paste(str1, str2, str3, str4, str5, sep = "<br/>" ))
      })
    
    
    
    
    
############### Figure resultat : Special Cases ##########
    
    output$result <- renderPlot({
      data <- figure()
      # print("Afficher les donnees : ")
      # print(data)
      
      n <- dim(data)[1]   # nb de donnees
      p <- dim(data)[2]   # nb d attributs ou dimension
      
      par(mfrow=c(1,2))
      
      fig <- data[,1:2]
      colnames(fig) <- c("x", "y")
      plot(fig, 
           xlim=c(0, 10), 
           ylim=c(0, 10),
           cex=1.3, pch=3, lwd=2)
      title(main=paste("Three Cases A, B and C ( n =", n, ") "))
      labels <- c("C","A","B")
      j <- 1
      for(i in 48:50){
        text(data[i,1], data[i,2], label=labels[j], cex= 1.5, pos=4)
        j <- j+1
        points(fig[i,1], fig[i,2], lty="solid", pch=3, cex=1.3, lwd=3, col="red")
      }
      draw.ellipse(5, 5, a=2, b=2, lty="dashed")
      draw.ellipse(5, 5, a=4, b=4, lty="dashed")
      
      lesObservations <- obs(data)
      # print("Afficher les observations : ")
      # print(lesObservations)
      
      top <- letop()
      debut_alpha <- 2
      fin_alpha <- n
      alpha <- detect_alpha(data, lesObservations, top, debut_alpha, fin_alpha)
      # print(" alpha = ")
      # print(alpha)
      # alpha <- 4
      isolAttribut <- affiche_rareness(data, lesObservations, alpha,  alpha)
      classement <- rank(-apply(isolAttribut, 1, max), ties.method = "first")
      letop <- which(classement<=top)
      
      plot(fig, 
           xlim=c(0, 10), 
           ylim=c(0, 10),
           cex=1.3, pch=3, lwd=2)
      title(main=paste("", top, "Special Cases ( with k = ", alpha, ")"))
      # title(sub=paste("Three special cases obtained using 5 nearest neighbours"))
      labels <- c("C","A","B")
      j <- 1
      for(i in 48:50){
        text(data[i,1], data[i,2], label=labels[j], cex= 1.5, pos=4)
        j <- j+1
        points(fig[i,1], fig[i,2], lty="solid", pch=3, cex=1.3, lwd=2, col="red")
      }
      
      for(i in letop){
        points(data[i,1], data[i,2], lty="solid", col="red", lwd=3, pch=1, cex=2)
      }
      coeffIsol <- apply(isolAttribut, 1, max)
      # print(" isolAttribut (rareness pour 2 voisinaage en ligne) :")
      # print(isolAttribut)
      # print(coeffIsol)
      for(i in letop){
        text(data[i,1], data[i,2], label=round(coeffIsol[i],4), cex= 1.5, pos=1)
      }
    })
    
    
      
})