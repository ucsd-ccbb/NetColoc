import concurrent.futures
import multiprocessing

def main():
  hiList = multiprocessing.RawArray('i', 10)
  print('start')
  with concurrent.futures.ProcessPoolExecutor() as executor:
    results = executor.map(simple, [[1,hiList], [3,hiList], [5,hiList]])
  print(hiList)
  for result in results:
    print(result)
  print(hiList)

def simple(args):
  return args


if __name__ == "__main__":
    main()

