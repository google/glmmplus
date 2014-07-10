library(RUnit)
library(glmmplus)

test.suite <- defineTestSuite("example",
                              dirs = getwd(),
                              testFileRegexp = '.*test\\.R$',
                              testFuncRegexp = 'Test.*')

test.result <- runTestSuite(test.suite)

printTextProtocol(test.result)
