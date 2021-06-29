!!Macro for assertions
#define ASSERTC(OBJ1,OPR,EXPECT) if(.not. (OBJ1 OPR EXPECT)) write(6,*) "ERROR"
#define ASSERTC2(OBJ1,OPR,EXPECT) if(.not. (OBJ1 OPR EXPECT)) write(6,*) "ERROR"
#define ASSERTCP(OBJ1,OPR,EXPECT) 
#define ASSERTI(bool,msg) if(.not. (bool)) write(6,*) msg

