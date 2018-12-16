#include<fstream>
#include<iostream>
#include<stdio.h>
#include<WinSock2.h>
#include<stdlib.h>
#include <time.h> 

using namespace std ;

#pragma comment (lib,"Ws2_32.lib")
 
int main()
{ 
	while (1)
	{
		WSADATA wsadata ;
		SOCKET  listening ;
		SOCKET  NEW ;
		SOCKADDR_IN  ServerAddr ;
		SOCKADDR_IN  KHD ;
		int KHD_len ;
		int port = 5150 ;
		int ret   ;
		//  初始化套接字
		WSAStartup(MAKEWORD(2,2),&wsadata) ;
		//  加载套接字
		listening  = socket(AF_INET,SOCK_STREAM,IPPROTO_TCP);
		
		ServerAddr.sin_family = AF_INET ;
		ServerAddr.sin_port = htons(port);
		ServerAddr.sin_addr.S_un.S_addr = htonl(INADDR_ANY);
		//	绑定套接字
		bind(listening,(SOCKADDR*)&ServerAddr,sizeof(ServerAddr));
		//	监听套接字
		listen(listening,50);
		//  连接
		KHD_len  = sizeof(SOCKADDR);
		NEW  = accept(listening,(SOCKADDR*)&KHD,&KHD_len);
		//  数据接收
		char ip[256] ;
		recv(NEW,ip,256,0); 
		cout << ip<< ": 连接成功..." << endl;
		char questionId[256] ;
		ret = recv(NEW,questionId,256,0);
		cout << questionId <<endl ;
		
		// 时间
		time_t rawtime; 
		struct tm * timeinfo; 
		time ( &rawtime ); 
		timeinfo = localtime ( &rawtime ); 
		printf ("提交时间: %s", asctime(timeinfo) ); 

		//  清理工作 
		closesocket(listening);
		closesocket(NEW);
		WSACleanup();
	}

	return 0 ;
}