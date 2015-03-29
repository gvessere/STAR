#include <cstdlib>

void HandleSignals(int s);
void ConfigureSignalHandlers(void (*cleanup)());
void RemoveSignalHandlers();
