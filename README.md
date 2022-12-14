# Optical-Amplifier-Erbium-Doped-Fiber-Amplifier-EDFA-

 Vidéo de la simulation :

[![Optical-Amplifier-Erbium-Doped-Fiber-Amplifier-EDFA](https://img.youtube.com/vi/HNaS2a-iUX0/0.jpg)](https://www.youtube.com/watch?v=HNaS2a-iUX0)

Les amplificateurs à fibre dopée erbium  sont des amplificateurs optiques qui utilisent une    fibre optique dopée comme support de gain pour amplifier un signal optique. Ils sont liés aux lasers à fibre. Le signal à amplifier et un laser de pompage sont multiplexés dans la fibre dopée, et le signal est amplifié par interaction avec les ions dopants comme le montre la figure suivante: ![imag1](https://user-images.githubusercontent.com/22806623/190588084-87236de8-47a1-4c61-a74b-cec65b7c673a.jpg)

Afin de faciliter la comparaison des différentes méthodes numériques et de faire le dimensionnement des amplificateurs à fibres dopées, nous sommes amenés à implémenter des interfaces graphiques permettant afin de prédire graphiquement et de façon ergonomique les performances d’un amplificateur EDFA avant même qu’il soit assemblé et caractérisé expérimentalement en laboratoire. A partir des paramètres qui seront accessibles au niveau de l’interface graphique, l’usager sera capable de simuler et si nécessaire comparer simultanément plusieurs configurations d’amplificateur par modification de : 

          * la puissance du laser de pompe
          * la longueur de la fibre dopée 
          * le type de pompage (pompage dans le même sens et sens contraire que la puissance du signal et le pompage dans les deux sens que la puissance du signal)
          * la nature de la fibre dopée  
          * le step des puissances(Ps, Pp )
          * le step de la longueur de la fibre etc.
       
Dans la suite de ce rapport nous détaillerons en premier point le principe



 de fonctionnement des  amplificateurs dopées à l'erbium et par la suite , nous détaillerons les étapes d'implémentation de l’interface graphique et en fin nous mettrons en applications l’interface graphique en cherchant les longueurs optimales des quatres différentes fibres : RadIX, Telcom, nLIGTH_Er80-8\128 et nLIGTH_Er110-4\128.
L'organigramme du GUI est illustré dans la figure suivante:
![gui_flow](https://user-images.githubusercontent.com/22806623/190591013-0bd1b078-5150-4ce7-9e3e-05a8898f2135.png)

La figure suivante represente les simulations obtenues:
![img5](https://user-images.githubusercontent.com/22806623/190592703-2172cd19-9c60-44be-bf0e-251e4cf10298.jpg)

![img7](https://user-images.githubusercontent.com/22806623/190592715-b04e9e01-8f0d-412f-9cc1-f67a7dad7e43.jpg)
