

```python
import pandas as pd
```


```python
blast_output_filename = "../data/results_99"
blast_output = pd.read_csv(blast_output_filename, sep="\t", header=None, names=[
    'query', 'subject', '%id', 'aln len', 'mismatches', 'gap opens', 'qstart', 'qend', 'sstart', 'send', 'eval', 
    'score'
], dtype={'query': 'object', 'subject': 'object', '%id': 'float32', 'aln len': 'int', 'mismatches': 'int', 
          'gap opens': 'int', 'qstart': 'int', 'qend': 'int', 'sstart': 'int', 'send': 'int', 
          'eval': 'float32', 'score': 'int'})
```


```python
blast_output.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>query</th>
      <th>subject</th>
      <th>%id</th>
      <th>aln len</th>
      <th>mismatches</th>
      <th>gap opens</th>
      <th>qstart</th>
      <th>qend</th>
      <th>sstart</th>
      <th>send</th>
      <th>eval</th>
      <th>score</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>NR_121333.1</td>
      <td>SRR5170882.38</td>
      <td>99.620003</td>
      <td>263</td>
      <td>1</td>
      <td>0</td>
      <td>265</td>
      <td>527</td>
      <td>1</td>
      <td>263</td>
      <td>0.0</td>
      <td>481</td>
    </tr>
    <tr>
      <th>1</th>
      <td>NR_121333.1</td>
      <td>SRR5170882.3</td>
      <td>98.859001</td>
      <td>263</td>
      <td>3</td>
      <td>0</td>
      <td>265</td>
      <td>527</td>
      <td>1</td>
      <td>263</td>
      <td>0.0</td>
      <td>470</td>
    </tr>
    <tr>
      <th>2</th>
      <td>NR_121333.1</td>
      <td>SRR5170882.49</td>
      <td>98.859001</td>
      <td>263</td>
      <td>2</td>
      <td>1</td>
      <td>265</td>
      <td>527</td>
      <td>1</td>
      <td>262</td>
      <td>0.0</td>
      <td>468</td>
    </tr>
    <tr>
      <th>3</th>
      <td>NR_121333.1</td>
      <td>SRR5170882.25</td>
      <td>98.859001</td>
      <td>263</td>
      <td>2</td>
      <td>1</td>
      <td>265</td>
      <td>527</td>
      <td>1</td>
      <td>262</td>
      <td>0.0</td>
      <td>468</td>
    </tr>
    <tr>
      <th>4</th>
      <td>NR_121332.1</td>
      <td>SRR5170882.38</td>
      <td>99.620003</td>
      <td>263</td>
      <td>1</td>
      <td>0</td>
      <td>265</td>
      <td>527</td>
      <td>1</td>
      <td>263</td>
      <td>0.0</td>
      <td>481</td>
    </tr>
  </tbody>
</table>
</div>




```python
from parsing_scripts.helper_functions import *
```


```python
blast_output2 = add_taxon_info(blast_output, '../../extra/fungi.ITS.fna')
blast_output2.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>query</th>
      <th>subject</th>
      <th>%id</th>
      <th>aln len</th>
      <th>mismatches</th>
      <th>gap opens</th>
      <th>qstart</th>
      <th>qend</th>
      <th>sstart</th>
      <th>send</th>
      <th>eval</th>
      <th>score</th>
      <th>taxid</th>
      <th>scientific_name</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>NR_121333.1</td>
      <td>SRR5170882.38</td>
      <td>99.620003</td>
      <td>263</td>
      <td>1</td>
      <td>0</td>
      <td>265</td>
      <td>527</td>
      <td>1</td>
      <td>263</td>
      <td>0.0</td>
      <td>481</td>
      <td>470172</td>
      <td>Cladosporium ossifragi CBS 842.91</td>
    </tr>
    <tr>
      <th>1</th>
      <td>NR_121333.1</td>
      <td>SRR5170882.3</td>
      <td>98.859001</td>
      <td>263</td>
      <td>3</td>
      <td>0</td>
      <td>265</td>
      <td>527</td>
      <td>1</td>
      <td>263</td>
      <td>0.0</td>
      <td>470</td>
      <td>470172</td>
      <td>Cladosporium ossifragi CBS 842.91</td>
    </tr>
    <tr>
      <th>2</th>
      <td>NR_121333.1</td>
      <td>SRR5170882.49</td>
      <td>98.859001</td>
      <td>263</td>
      <td>2</td>
      <td>1</td>
      <td>265</td>
      <td>527</td>
      <td>1</td>
      <td>262</td>
      <td>0.0</td>
      <td>468</td>
      <td>470172</td>
      <td>Cladosporium ossifragi CBS 842.91</td>
    </tr>
    <tr>
      <th>3</th>
      <td>NR_121333.1</td>
      <td>SRR5170882.25</td>
      <td>98.859001</td>
      <td>263</td>
      <td>2</td>
      <td>1</td>
      <td>265</td>
      <td>527</td>
      <td>1</td>
      <td>262</td>
      <td>0.0</td>
      <td>468</td>
      <td>470172</td>
      <td>Cladosporium ossifragi CBS 842.91</td>
    </tr>
    <tr>
      <th>4</th>
      <td>NR_121332.1</td>
      <td>SRR5170882.38</td>
      <td>99.620003</td>
      <td>263</td>
      <td>1</td>
      <td>0</td>
      <td>265</td>
      <td>527</td>
      <td>1</td>
      <td>263</td>
      <td>0.0</td>
      <td>481</td>
      <td>470179</td>
      <td>Cladosporium antarcticum CBS 690.92</td>
    </tr>
  </tbody>
</table>
</div>




```python
blast_output['accession'] = blast_output.apply(lambda x: x['subject'].split('.')[0], axis=1)
blast_output['combined'] = blast_output[['query', 'accession']].apply(lambda x: ' '.join(x), axis=1)
blast_output2 = blast_output.groupby(["query", "accession"])["combined"].count().reset_index(name="count")
blast_output2.head(20)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>query</th>
      <th>accession</th>
      <th>count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>NR_073209.1</td>
      <td>DRR004164</td>
      <td>26</td>
    </tr>
    <tr>
      <th>1</th>
      <td>NR_073209.1</td>
      <td>DRR006809</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>NR_073209.1</td>
      <td>DRR006980</td>
      <td>2</td>
    </tr>
    <tr>
      <th>3</th>
      <td>NR_073209.1</td>
      <td>DRR007269</td>
      <td>6</td>
    </tr>
    <tr>
      <th>4</th>
      <td>NR_073209.1</td>
      <td>DRR010776</td>
      <td>1</td>
    </tr>
    <tr>
      <th>5</th>
      <td>NR_073209.1</td>
      <td>ERR1120559</td>
      <td>1</td>
    </tr>
    <tr>
      <th>6</th>
      <td>NR_073209.1</td>
      <td>ERR1120585</td>
      <td>1</td>
    </tr>
    <tr>
      <th>7</th>
      <td>NR_073209.1</td>
      <td>ERR1120614</td>
      <td>1</td>
    </tr>
    <tr>
      <th>8</th>
      <td>NR_073209.1</td>
      <td>ERR1120646</td>
      <td>1</td>
    </tr>
    <tr>
      <th>9</th>
      <td>NR_073209.1</td>
      <td>ERR1120667</td>
      <td>1</td>
    </tr>
    <tr>
      <th>10</th>
      <td>NR_073209.1</td>
      <td>ERR1456960</td>
      <td>1</td>
    </tr>
    <tr>
      <th>11</th>
      <td>NR_073209.1</td>
      <td>ERR1686646</td>
      <td>14</td>
    </tr>
    <tr>
      <th>12</th>
      <td>NR_073209.1</td>
      <td>ERR1701229</td>
      <td>1</td>
    </tr>
    <tr>
      <th>13</th>
      <td>NR_073209.1</td>
      <td>ERR1701242</td>
      <td>1</td>
    </tr>
    <tr>
      <th>14</th>
      <td>NR_073209.1</td>
      <td>ERR1742934</td>
      <td>3</td>
    </tr>
    <tr>
      <th>15</th>
      <td>NR_073209.1</td>
      <td>ERR1949666</td>
      <td>1</td>
    </tr>
    <tr>
      <th>16</th>
      <td>NR_073209.1</td>
      <td>ERR1949714</td>
      <td>1</td>
    </tr>
    <tr>
      <th>17</th>
      <td>NR_073209.1</td>
      <td>ERR1993823</td>
      <td>1</td>
    </tr>
    <tr>
      <th>18</th>
      <td>NR_073209.1</td>
      <td>ERR1993892</td>
      <td>1</td>
    </tr>
    <tr>
      <th>19</th>
      <td>NR_073209.1</td>
      <td>ERR1993925</td>
      <td>2</td>
    </tr>
  </tbody>
</table>
</div>




```python
blast_output3 = blast_output3.iloc[1:]
blast_output3.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>query</th>
      <th>accession</th>
      <th>count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1</th>
      <td>NR_073209.1</td>
      <td>DRR004164</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>NR_073209.1</td>
      <td>DRR006809</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>NR_073209.1</td>
      <td>DRR006980</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>NR_073209.1</td>
      <td>DRR007269</td>
      <td>1</td>
    </tr>
    <tr>
      <th>5</th>
      <td>NR_073209.1</td>
      <td>DRR010776</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>




```python
pd.options.mode.chained_assignment = None
```


```python
blast_output2['combined'] = blast_output2[['query', 'accession']].apply(lambda x: ' '.join(x), axis=1)
blast_output3 = blast_output2.groupby(["query", "accession"])["combined"].count().reset_index(name="count")
blast_output3.head(20)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>query</th>
      <th>accession</th>
      <th>count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>NR_073209.1</td>
      <td>DRR004164</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>NR_073209.1</td>
      <td>DRR006809</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>NR_073209.1</td>
      <td>DRR006980</td>
      <td>1</td>
    </tr>
    <tr>
      <th>3</th>
      <td>NR_073209.1</td>
      <td>DRR007269</td>
      <td>1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>NR_073209.1</td>
      <td>DRR010776</td>
      <td>1</td>
    </tr>
    <tr>
      <th>5</th>
      <td>NR_073209.1</td>
      <td>ERR1120559</td>
      <td>1</td>
    </tr>
    <tr>
      <th>6</th>
      <td>NR_073209.1</td>
      <td>ERR1120585</td>
      <td>1</td>
    </tr>
    <tr>
      <th>7</th>
      <td>NR_073209.1</td>
      <td>ERR1120614</td>
      <td>1</td>
    </tr>
    <tr>
      <th>8</th>
      <td>NR_073209.1</td>
      <td>ERR1120646</td>
      <td>1</td>
    </tr>
    <tr>
      <th>9</th>
      <td>NR_073209.1</td>
      <td>ERR1120667</td>
      <td>1</td>
    </tr>
    <tr>
      <th>10</th>
      <td>NR_073209.1</td>
      <td>ERR1456960</td>
      <td>1</td>
    </tr>
    <tr>
      <th>11</th>
      <td>NR_073209.1</td>
      <td>ERR1686646</td>
      <td>1</td>
    </tr>
    <tr>
      <th>12</th>
      <td>NR_073209.1</td>
      <td>ERR1701229</td>
      <td>1</td>
    </tr>
    <tr>
      <th>13</th>
      <td>NR_073209.1</td>
      <td>ERR1701242</td>
      <td>1</td>
    </tr>
    <tr>
      <th>14</th>
      <td>NR_073209.1</td>
      <td>ERR1742934</td>
      <td>1</td>
    </tr>
    <tr>
      <th>15</th>
      <td>NR_073209.1</td>
      <td>ERR1949666</td>
      <td>1</td>
    </tr>
    <tr>
      <th>16</th>
      <td>NR_073209.1</td>
      <td>ERR1949714</td>
      <td>1</td>
    </tr>
    <tr>
      <th>17</th>
      <td>NR_073209.1</td>
      <td>ERR1993823</td>
      <td>1</td>
    </tr>
    <tr>
      <th>18</th>
      <td>NR_073209.1</td>
      <td>ERR1993892</td>
      <td>1</td>
    </tr>
    <tr>
      <th>19</th>
      <td>NR_073209.1</td>
      <td>ERR1993925</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>




```python
blast_output3 = blast_output2.groupby(["query", "accession"])["combined"].count().reset_index(name="count")
blast_output3.head(20)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>query</th>
      <th>accession</th>
      <th>count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>NR_073208.1</td>
      <td>DRR004164</td>
      <td>29</td>
    </tr>
    <tr>
      <th>1</th>
      <td>NR_073208.1</td>
      <td>DRR004352</td>
      <td>1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>NR_073208.1</td>
      <td>DRR004922</td>
      <td>2</td>
    </tr>
    <tr>
      <th>3</th>
      <td>NR_073208.1</td>
      <td>DRR004964</td>
      <td>2</td>
    </tr>
    <tr>
      <th>4</th>
      <td>NR_073208.1</td>
      <td>DRR010490</td>
      <td>1</td>
    </tr>
    <tr>
      <th>5</th>
      <td>NR_073208.1</td>
      <td>DRR020111</td>
      <td>1</td>
    </tr>
    <tr>
      <th>6</th>
      <td>NR_073208.1</td>
      <td>ERR1120559</td>
      <td>1</td>
    </tr>
    <tr>
      <th>7</th>
      <td>NR_073208.1</td>
      <td>ERR1120585</td>
      <td>1</td>
    </tr>
    <tr>
      <th>8</th>
      <td>NR_073208.1</td>
      <td>ERR1120598</td>
      <td>1</td>
    </tr>
    <tr>
      <th>9</th>
      <td>NR_073208.1</td>
      <td>ERR1120667</td>
      <td>1</td>
    </tr>
    <tr>
      <th>10</th>
      <td>NR_073208.1</td>
      <td>ERR1173801</td>
      <td>5</td>
    </tr>
    <tr>
      <th>11</th>
      <td>NR_073208.1</td>
      <td>ERR1279121</td>
      <td>2</td>
    </tr>
    <tr>
      <th>12</th>
      <td>NR_073208.1</td>
      <td>ERR1456960</td>
      <td>210</td>
    </tr>
    <tr>
      <th>13</th>
      <td>NR_073208.1</td>
      <td>ERR1993765</td>
      <td>1</td>
    </tr>
    <tr>
      <th>14</th>
      <td>NR_073208.1</td>
      <td>ERR1993823</td>
      <td>1</td>
    </tr>
    <tr>
      <th>15</th>
      <td>NR_073208.1</td>
      <td>ERR1993892</td>
      <td>1</td>
    </tr>
    <tr>
      <th>16</th>
      <td>NR_073208.1</td>
      <td>ERR1993925</td>
      <td>2</td>
    </tr>
    <tr>
      <th>17</th>
      <td>NR_073208.1</td>
      <td>ERR2615922</td>
      <td>1</td>
    </tr>
    <tr>
      <th>18</th>
      <td>NR_073208.1</td>
      <td>ERR677491</td>
      <td>8</td>
    </tr>
    <tr>
      <th>19</th>
      <td>NR_073208.1</td>
      <td>ERR677497</td>
      <td>22</td>
    </tr>
  </tbody>
</table>
</div>




```python
blast_output3.rename(columns={'count':'hits', 'accession': 'subject_accession'}, inplace=True)
blast_output4 = add_taxon_info(blast_output3, '../../extra/fungi.ITS.fna')
blast_output4.head(20)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>query</th>
      <th>subject_accession</th>
      <th>hits</th>
      <th>taxid</th>
      <th>scientific_name</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>NR_073209.1</td>
      <td>DRR004164</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
    </tr>
    <tr>
      <th>1</th>
      <td>NR_073209.1</td>
      <td>DRR006809</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
    </tr>
    <tr>
      <th>2</th>
      <td>NR_073209.1</td>
      <td>DRR006980</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
    </tr>
    <tr>
      <th>3</th>
      <td>NR_073209.1</td>
      <td>DRR007269</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
    </tr>
    <tr>
      <th>4</th>
      <td>NR_073209.1</td>
      <td>DRR010776</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
    </tr>
    <tr>
      <th>5</th>
      <td>NR_073209.1</td>
      <td>ERR1120559</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
    </tr>
    <tr>
      <th>6</th>
      <td>NR_073209.1</td>
      <td>ERR1120585</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
    </tr>
    <tr>
      <th>7</th>
      <td>NR_073209.1</td>
      <td>ERR1120614</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
    </tr>
    <tr>
      <th>8</th>
      <td>NR_073209.1</td>
      <td>ERR1120646</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
    </tr>
    <tr>
      <th>9</th>
      <td>NR_073209.1</td>
      <td>ERR1120667</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
    </tr>
    <tr>
      <th>10</th>
      <td>NR_073209.1</td>
      <td>ERR1456960</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
    </tr>
    <tr>
      <th>11</th>
      <td>NR_073209.1</td>
      <td>ERR1686646</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
    </tr>
    <tr>
      <th>12</th>
      <td>NR_073209.1</td>
      <td>ERR1701229</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
    </tr>
    <tr>
      <th>13</th>
      <td>NR_073209.1</td>
      <td>ERR1701242</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
    </tr>
    <tr>
      <th>14</th>
      <td>NR_073209.1</td>
      <td>ERR1742934</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
    </tr>
    <tr>
      <th>15</th>
      <td>NR_073209.1</td>
      <td>ERR1949666</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
    </tr>
    <tr>
      <th>16</th>
      <td>NR_073209.1</td>
      <td>ERR1949714</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
    </tr>
    <tr>
      <th>17</th>
      <td>NR_073209.1</td>
      <td>ERR1993823</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
    </tr>
    <tr>
      <th>18</th>
      <td>NR_073209.1</td>
      <td>ERR1993892</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
    </tr>
    <tr>
      <th>19</th>
      <td>NR_073209.1</td>
      <td>ERR1993925</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
    </tr>
  </tbody>
</table>
</div>




```python
db = SRAdb('SRAmetadb.sqlite')
```


```python
output = add_biosample_title(blast_output4.head(500), db)
```


```python
def add_biosample_title(df, db):
    """The 'db' input must be an SRAdb object (from the pysradb library) connected to a copy of the SRAmetadb.sqlite database."""
    df['biosample_title'] = df.apply(lambda x: get_experiment_title(x["subject_accession"], db).get("experiment_title")[0], axis=1)
    return df

def add_biosample_attribute(df, db):
    """The 'db' input must be an SRAdb object (from the pysradb library) connected to a copy of the SRAmetadb.sqlite database."""
    df['biosample_attr'] = df.apply(lambda x: get_experiment_title(x["subject_accession"], db).get("sample_attribute")[0], axis=1)
    return df
```


```python
output
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>query</th>
      <th>subject_accession</th>
      <th>hits</th>
      <th>taxid</th>
      <th>scientific_name</th>
      <th>biosample_title</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>NR_073209.1</td>
      <td>DRR004164</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>YamamotoYoshidayama__ACTGAGTC__ITS</td>
    </tr>
    <tr>
      <th>1</th>
      <td>NR_073209.1</td>
      <td>DRR006809</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>HVL5RMY01__TOMA_048__ITS3</td>
    </tr>
    <tr>
      <th>2</th>
      <td>NR_073209.1</td>
      <td>DRR006980</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>HVL5RMY01__TOMA_222__ITS3</td>
    </tr>
    <tr>
      <th>3</th>
      <td>NR_073209.1</td>
      <td>DRR007269</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>H33YU2J01__TOMA_520__ITS3</td>
    </tr>
    <tr>
      <th>4</th>
      <td>NR_073209.1</td>
      <td>DRR010776</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>HUMWT9A01__Yaku_0428__ITS3</td>
    </tr>
    <tr>
      <th>5</th>
      <td>NR_073209.1</td>
      <td>ERR1120559</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>454 GS FLX Titanium sequencing</td>
    </tr>
    <tr>
      <th>6</th>
      <td>NR_073209.1</td>
      <td>ERR1120585</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>454 GS FLX+ sequencing</td>
    </tr>
    <tr>
      <th>7</th>
      <td>NR_073209.1</td>
      <td>ERR1120614</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>454 GS FLX Titanium sequencing</td>
    </tr>
    <tr>
      <th>8</th>
      <td>NR_073209.1</td>
      <td>ERR1120646</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>454 GS FLX Titanium sequencing</td>
    </tr>
    <tr>
      <th>9</th>
      <td>NR_073209.1</td>
      <td>ERR1120667</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>454 GS FLX Titanium sequencing</td>
    </tr>
    <tr>
      <th>10</th>
      <td>NR_073209.1</td>
      <td>ERR1456960</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>454 GS FLX sequencing</td>
    </tr>
    <tr>
      <th>11</th>
      <td>NR_073209.1</td>
      <td>ERR1686646</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>None</td>
    </tr>
    <tr>
      <th>12</th>
      <td>NR_073209.1</td>
      <td>ERR1701229</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>Illumina MiSeq sequencing</td>
    </tr>
    <tr>
      <th>13</th>
      <td>NR_073209.1</td>
      <td>ERR1701242</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>Illumina MiSeq sequencing</td>
    </tr>
    <tr>
      <th>14</th>
      <td>NR_073209.1</td>
      <td>ERR1742934</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>Illumina MiSeq paired end sequencing</td>
    </tr>
    <tr>
      <th>15</th>
      <td>NR_073209.1</td>
      <td>ERR1949666</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>None</td>
    </tr>
    <tr>
      <th>16</th>
      <td>NR_073209.1</td>
      <td>ERR1949714</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>None</td>
    </tr>
    <tr>
      <th>17</th>
      <td>NR_073209.1</td>
      <td>ERR1993823</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>None</td>
    </tr>
    <tr>
      <th>18</th>
      <td>NR_073209.1</td>
      <td>ERR1993892</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>None</td>
    </tr>
    <tr>
      <th>19</th>
      <td>NR_073209.1</td>
      <td>ERR1993925</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>None</td>
    </tr>
    <tr>
      <th>20</th>
      <td>NR_073209.1</td>
      <td>ERR2107844</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>None</td>
    </tr>
    <tr>
      <th>21</th>
      <td>NR_073209.1</td>
      <td>ERR2107849</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>None</td>
    </tr>
    <tr>
      <th>22</th>
      <td>NR_073209.1</td>
      <td>ERR2107853</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>None</td>
    </tr>
    <tr>
      <th>23</th>
      <td>NR_073209.1</td>
      <td>ERR2107854</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>None</td>
    </tr>
    <tr>
      <th>24</th>
      <td>NR_073209.1</td>
      <td>ERR2107855</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>None</td>
    </tr>
    <tr>
      <th>25</th>
      <td>NR_073209.1</td>
      <td>ERR2107859</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>None</td>
    </tr>
    <tr>
      <th>26</th>
      <td>NR_073209.1</td>
      <td>ERR2107860</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>None</td>
    </tr>
    <tr>
      <th>27</th>
      <td>NR_073209.1</td>
      <td>ERR2107869</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>None</td>
    </tr>
    <tr>
      <th>28</th>
      <td>NR_073209.1</td>
      <td>ERR2107873</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>None</td>
    </tr>
    <tr>
      <th>29</th>
      <td>NR_073209.1</td>
      <td>ERR2107876</td>
      <td>1</td>
      <td>105984</td>
      <td>Apiotrichum porosum CBS 2040</td>
      <td>None</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>470</th>
      <td>NR_073212.1</td>
      <td>SRR3945903</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>ITS targeted metagenome: fungal community from...</td>
    </tr>
    <tr>
      <th>471</th>
      <td>NR_073212.1</td>
      <td>SRR3945913</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>ITS targeted metagenome: fungal community from...</td>
    </tr>
    <tr>
      <th>472</th>
      <td>NR_073212.1</td>
      <td>SRR3945916</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>ITS targeted metagenome: fungal community from...</td>
    </tr>
    <tr>
      <th>473</th>
      <td>NR_073212.1</td>
      <td>SRR3945918</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>ITS targeted metagenome: fungal community from...</td>
    </tr>
    <tr>
      <th>474</th>
      <td>NR_073212.1</td>
      <td>SRR3945919</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>ITS targeted metagenome: fungal community from...</td>
    </tr>
    <tr>
      <th>475</th>
      <td>NR_073212.1</td>
      <td>SRR3945946</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>ITS targeted metagenome: fungal community from...</td>
    </tr>
    <tr>
      <th>476</th>
      <td>NR_073212.1</td>
      <td>SRR3945951</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>ITS targeted metagenome: fungal community from...</td>
    </tr>
    <tr>
      <th>477</th>
      <td>NR_073212.1</td>
      <td>SRR3945959</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>ITS targeted metagenome: fungal community from...</td>
    </tr>
    <tr>
      <th>478</th>
      <td>NR_073212.1</td>
      <td>SRR4308344</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>rRNA ITS2 amplicons with fungal specific primers</td>
    </tr>
    <tr>
      <th>479</th>
      <td>NR_073212.1</td>
      <td>SRR4308396</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>rRNA ITS2 amplicons with fungal specific primers</td>
    </tr>
    <tr>
      <th>480</th>
      <td>NR_073212.1</td>
      <td>SRR4308404</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>rRNA ITS2 amplicons with fungal specific primers</td>
    </tr>
    <tr>
      <th>481</th>
      <td>NR_073212.1</td>
      <td>SRR4308405</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>rRNA ITS2 amplicons with fungal specific primers</td>
    </tr>
    <tr>
      <th>482</th>
      <td>NR_073212.1</td>
      <td>SRR4308410</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>rRNA ITS2 amplicons with fungal specific primers</td>
    </tr>
    <tr>
      <th>483</th>
      <td>NR_073212.1</td>
      <td>SRR5000090</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>None</td>
    </tr>
    <tr>
      <th>484</th>
      <td>NR_073212.1</td>
      <td>SRR5000130</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>None</td>
    </tr>
    <tr>
      <th>485</th>
      <td>NR_073212.1</td>
      <td>SRR5361826</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>None</td>
    </tr>
    <tr>
      <th>486</th>
      <td>NR_073212.1</td>
      <td>SRR5361829</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>None</td>
    </tr>
    <tr>
      <th>487</th>
      <td>NR_073212.1</td>
      <td>SRR5361841</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>None</td>
    </tr>
    <tr>
      <th>488</th>
      <td>NR_073212.1</td>
      <td>SRR5361844</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>None</td>
    </tr>
    <tr>
      <th>489</th>
      <td>NR_073212.1</td>
      <td>SRR5361850</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>None</td>
    </tr>
    <tr>
      <th>490</th>
      <td>NR_073212.1</td>
      <td>SRR5361906</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>None</td>
    </tr>
    <tr>
      <th>491</th>
      <td>NR_073212.1</td>
      <td>SRR5361919</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>None</td>
    </tr>
    <tr>
      <th>492</th>
      <td>NR_073212.1</td>
      <td>SRR5364708</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>None</td>
    </tr>
    <tr>
      <th>493</th>
      <td>NR_073212.1</td>
      <td>SRR5364739</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>None</td>
    </tr>
    <tr>
      <th>494</th>
      <td>NR_073212.1</td>
      <td>SRR5364786</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>None</td>
    </tr>
    <tr>
      <th>495</th>
      <td>NR_073212.1</td>
      <td>SRR5383979</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>None</td>
    </tr>
    <tr>
      <th>496</th>
      <td>NR_073212.1</td>
      <td>SRR5573690</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>None</td>
    </tr>
    <tr>
      <th>497</th>
      <td>NR_073212.1</td>
      <td>SRR5573693</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>None</td>
    </tr>
    <tr>
      <th>498</th>
      <td>NR_073212.1</td>
      <td>SRR5573700</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>None</td>
    </tr>
    <tr>
      <th>499</th>
      <td>NR_073212.1</td>
      <td>SRR5573701</td>
      <td>1</td>
      <td>40409</td>
      <td>Solicoccozyma terreus CBS 1895</td>
      <td>None</td>
    </tr>
  </tbody>
</table>
<p>500 rows × 6 columns</p>
</div>


