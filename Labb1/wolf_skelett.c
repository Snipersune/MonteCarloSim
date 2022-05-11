int update(Par *par, int *spin){
    
    // En update-körning består av att bygga ETT 'cluster' av noder genom att använda Wolff cluster 
    // algoritmen, och flippa spinen av noderna man lägger till.

    accept <- 0 // För att kolla hur många noder som lagts till i clusteret.


    /* Wolf cluster algoritmen */
    start_index <- random int; // Slumpa startnod att börja bygga clusteret ifrån.
    start_state <- spin[start_index] // Spara startnodens spinstate. Detta är det state resten
                                     // av noderna i algoritmen kommer utvärderas mot.
    
    spin[start_index] *= -1 // Flippa state av startnod. 

    // Algoritmen fungerar sedan genom att gå igenom grannarna till noderna i clusteret, och 
    // kolla ifall de har samma spin som 'start_state'. Har de inte det ignorerar man dessa. 
    // Har de samma state som 'start_state' ska de läggas till i clusteret med sannolikhet 'prob'. 
    // Dessa nytillagda noders grannar måste då också kollas, och så fortsätter algoritmen tills dess
    // att inga nya noder läggs till.
    prob <- 1.0 - exp(-2.0/par->t)
    
    // Använd kö för att lägga till noder som ska kollas
    enqueue <- start_index // Lägg in start_index som första nod att kolla.


    // Här sker själva algoritmen som bygger upp clusteret.
    while queue not empty{
        accept++ // Kön va inte tom, ytterligare nod finns i clusteret

        node_index <- dequeue // Hämta index till nod att kolla och ta bort index från kön.
        nbors <- neighbours[node_index] // Hämta på något sätt index till nodens grannar.

        for nbor in nbors { // Gå igenom grannar
            if spin[nbor] == start_state { // Kolla ifall granne har rätt state
                if rand < prob { // Lägg till med sannolikhet 'prob'

                    spin[nbor] *= -1 // Noden kom med i clusteret, flippa spin.
                    enqueue <- nbor  // Lägg till grannens index i kön, eftersom noden kom med i
                                     // clusteret måste även denna nods grannar kollas.
                }
            }
        }
    }
    return accept
}

// Detta är i princip allt. Kön kommer fyllas på med de grannar som läggs till i 
// clusteret, och när deras grannar i sin tur kollas tas de bort från kön. När kön 
// sedan är tom är clusteret färdigbyggt, och alla noder där i kommer att ha flippat spin.

// Viktigt att flippa spinen på en nod så fort man lägger till den till clusteret, 
// annars om man har otur kan den läggas till fler gånger då den kan vara granne till
// andra noder man kommer kolla igenom också. 