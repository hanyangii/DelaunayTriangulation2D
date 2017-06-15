1. 실행 방법 
'Point.txt' 파일에 x, y좌표를 탭으로 구분하여 입력한 뒤 
DCDelaunay.c 소스를 실행하면 됩니다. 
중복된 점의 좌표를 살짝 옆으로 옮기는 방법을 이용하기 때문에 
'newPoints.txt'에 새로운 점들의 리스트가 생깁니다. 
헤더파일의 display 변수를 true로 지정후 실행하면 'Cout.txt' 파일을 통해 
Triangulation과정을 확인 하실 수 있습니다. 
결과는 'Link.txt'파일로 edge들의 리스트가 출력됩니다. 

2. 테스트 방법
'Delaunay_test.m' 파일의  'Link.txt'파일과 'newPoints.txt'파일의 경로를 지정해 준 뒤
실행하면 matlab의 결과와 본 소스에서의 결과를 비교하실 수 있습니다. 