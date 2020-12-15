#include "Louie3DEngine.h"
int main()
{
	Louie3DEngine demo;
	if (demo.Construct(800, 450, 2, 2))
		demo.Start();
    return 0;
}
