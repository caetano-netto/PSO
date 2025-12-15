Passos para Executar o Programa no Windows (CMD / PowerShell)

1️⃣ Abrir o terminal

Abra o Prompt de Comando (CMD) ou Windows PowerShell.

2️⃣ Acessar a pasta do projeto

Navegue até a pasta onde os arquivos demo.c e pso.c estão salvos.

cd caminho\para\pso-master


3️⃣ Compilar o programa

Compile o código com o GCC:

gcc demo.c pso.c -O2 -lm -o demo


4️⃣ Ajustar o terminal para UTF-8 (opcional, recomendado)

Para que os caracteres gráficos do menu apareçam corretamente, digite:

chcp 65001