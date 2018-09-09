---
layout: article
date: 2017-09-22 09:42
title: "ポインタのキャストがわかりにくい"
tags: [c, pointer, cast]
comments: true
published: true
---

この前，callback関数についてちょっと調べようと思って[Wikipediaのページ](https://ja.wikipedia.org/wiki/%E3%82%B3%E3%83%BC%E3%83%AB%E3%83%90%E3%83%83%E3%82%AF_(%E6%83%85%E5%A0%B1%E5%B7%A5%E5%AD%A6))のコードを見ていたが，ポインタ変数への代入部分でハマったのでメモする．

コールバックの例のコード (C言語):

```c
#include <stdio.h>

/* Library code */
int traverseWith(int array[], size_t length,
                 int (*callback)(int index, int item, void *param),
                 void *param)
{
    int exitCode = 0;
    for (int i = 0; i < length; i++) {
        exitCode = callback(i, array[i], param);
        if (exitCode) {
            break;
        }
    }
    return exitCode;
}

/* Application code */
int search (int index, int item, void *param)
{
    if (item > 5) {
        *(int *)param = index;
        return 1;
    } else {
        return 0;
    }
}

int main(int argc, char const* argv[])
{
  int index;
  int found;
  int array[6] = {1,2,3,4,5,6};
  int length = sizeof(array) / sizeof(array[0]);
  found = traverseWith(array, length, search, &index);
  if (found) {
      printf("Item %d\n", index);
  } else {
      printf("Not found\n");
  }

  return 0;
}
```

問題の部分はアプリケーションコードの以下の部分．

```c
    if (item > 5) {
        *(int *)param = index;
        return 1;
    }
```

`*(int *)param = index;`ってなんだ．．．？`param`は`(void *)`型の変数だから`*`をつけてポインタの参照先に`index`の値を代入しているんだろうが，`(int *)`がわからん．．．ポインタのポインタ？でもparamの値は`&index`で`int`の番地だし．．．

しばらく考えて分からなかったので，`*(int *)param = index;`の左辺を適当に変えてコンパイルしてみることにした．

結果，`(int *)`は元が`(void *)`である変数`param`のキャストで，頭の`*`がポインタの参照先を示す記号であることがわかった．  
よって`index`が代入されるメモリ番地は`main`内の`index`変数用に確保されている領域の番地なので，`main`内の`index`変数に`search`関数の結果である配列のインデックスが格納される．
