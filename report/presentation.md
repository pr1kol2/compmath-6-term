---
marp: true
theme: default
paginate: true
math: katex
size: 16:9
style: |
  section {
    font-size: 22px;
  }
  h1 {
    font-size: 36px;
    color: #1a237e;
  }
  h2 {
    font-size: 28px;
    color: #354199;
  }
  table {
    font-size: 18px;
  }
  .columns {
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 1em;
  }
  blockquote {
    border-left: 4px solid #1a237e;
    padding-left: 1em;
    font-style: italic;
    color: #555;
  }
  footer {
    font-size: 14px;
    color: #888;
  }
---

<!-- _class: lead -->
<!-- _paginate: false -->

# Применения интерполяции с регуляризацией для регрессии зашумлённых данных

## Кирилл Гринько Б05-332

---

# Постановка задачи

Пусть $x_1,\dots,x_m \in [a,b]$ — точки наблюдений, $y_i \in \mathbb{R}$ — измерения:

$$y_i = f(x_i) + \varepsilon_i, \quad i = 1, \dots, m$$

- $f \in C[a,b]$ — неизвестная гладкая функция,
- $\varepsilon_i \in \mathbb{R}$ — шум, $\mathbb{E}[\varepsilon_i]=0$, $\mathbb{D}[\varepsilon_i] = \sigma^2$

**Задача:** по наблюдениям $(x_i, y_i)$ построить оценку $\hat{f} \colon [a,b] \to \mathbb{R}$

**Проблема:** точная интерполяция $\hat{f}(x_i) = y_i$ проходит через зашумлённые данные — **переобучение** (overfitting).

---

# Константа Лебега

Работаем в пространстве $\bigl(C[a,b],\,\|\cdot\|_С\bigr)$, где $\|f\|_С = \max_{x\in[a,b]}|f(x)|$.

Пусть $x_0,\dots,x_n \in [a,b]$ — узлы интерполяции, $\Pi_n$ — пространство полиномов степени $\le n$.

Оператор интерполяции (проектор) $P_n \colon C[a,b] \to \Pi_n$ задаётся условиями $[P_n f](x_j) = f(x_j)$, $j=0,\dots,n$.

**Определение.** Константа Лебега — норма оператора $P_n$:
$$\Lambda_n \;=\; \|P_n\|_{C[a,b]\to C[a,b]} \;=\; \|P_n\| \;=\; \sup_{\substack{f \in C[a,b] \\ \|f\|_С=1}}\|P_n f\|_С$$

**Вычисление.** Из формы Лагранжа $[P_n f](x) = \sum_{j=0}^n f(x_j)\,l_j(x)$, где $l_j \in \Pi_n$ — базисные полиномы Лагранжа:
$$|[P_n f](x)| \;\le\; \|f\|_С \sum_{j=0}^n |l_j(x)| \;\Rightarrow\; \|P_n f\|_С \;\le\; \|f\|_С\cdot\max_{x\in[a,b]}\sum_{j=0}^n|l_j(x)|$$

---

# Константа Лебега

Равенство достигается на $f^*(x_j) = \operatorname{sign} l_j(x^*)$, $x^* = \operatorname{argmax}_{x}\sum_j|l_j(x)|$, поэтому:
$$\Lambda_n \;=\; \max_{x\in[a,b]}\sum_{j=0}^n |l_j(x)|$$

**Лемма Лебега.** Пусть $E_n^*(f) = \inf_{q\in\Pi_n}\|f-q\|_С$ — наилучшее приближение $f$ в $\bigl(\Pi_n, \|\cdot\|_С\bigr)$. Тогда:
$$\|f - P_n f\|_С \;\le\; (1 + \Lambda_n)\cdot E_n^*(f)$$

---

# Устойчивость к шуму в узловых значениях

Пусть узловые значения заданы с ошибкой:
$$y_j = f_j + \delta f_j, \qquad |\delta f_j| \le \varepsilon_j \le \delta, \quad j = 0,\dots,n$$

**Оценка устойчивости.** Из линейности $P_n$ и определения $\Lambda_n$:
$$\|P_ny - P_n f\|_C \;=\; \|P_n(\delta f)\|_C \;\le\; \|P_n\| \cdot \|\delta f\|_C \;\le\; \Lambda_n \cdot \delta$$

Равенство все также достигается на $f^*(x_j) = \operatorname{sign} l_j(x^*)$.

**Полная ошибка**:
$$\|f - P_ny\|_C \;\le\; \underbrace{(1+\Lambda_n)\,E_n^*(f)}_{\text{ошибка приближения}} \;+\; \underbrace{\Lambda_n\cdot\delta}_{\text{ошибка из-за шума}}$$

---

# Рост константы Лебега для различных сеток

Рассматриваем $P_n$ в $\bigl(C[-1,1],\,\|\cdot\|_\infty\bigr)$, сетка из $n+1$ узлов на $[-1,1]$.

| Сетка | Оценка |
|---|---|
| **Любая** (Теорема Бернштейна–Фабера) | $\Lambda_n \;\ge\; c \, ln \, n$ |
| **Равномерная** | $\Lambda_n \ge \dfrac{\theta \; (1 - \theta) \; 2^{n}}{n ^ 2}, \quad  0 < \theta < 1$ |
| **Чебышёвская** | $\Lambda_n = O(ln \, n)$ |

---

# Принцип регуляризации

**Общий принцип (Тихонов, 1963).** Вместо некорректной задачи точной интерполяции ищем $\hat{f}$ из пространства $\mathcal{F}$, минимизирующую функционал в $\mathbb{R}$:

$$\boxed{\mathcal{J}(\hat{f}) = \underbrace{\mathcal{D}(\hat{f},\, y)}_{\text{близость к данным}} + \lambda \cdot \underbrace{\Omega(\hat{f})}_{\text{штраф за сложность}}, \qquad \lambda > 0}$$

---

# Принцип регуляризации

**Стандартный выбор** $\mathcal{D}$: среднеквадратичная невязка $\mathcal{D}(\hat{f}, y) = \frac{1}{m}\sum_{i=1}^m (\hat{f}(x_i) - y_i)^2$.

**Выбор $\Omega$ определяется природой «нежелательного» поведения:**

- Если $\mathcal{F} = \Pi_n$ и решение параметризовано вектором коэффициентов $c \in \mathbb{R}^{n+1}$, то осцилляции полинома порождаются большими $\|c\|_2$. Естественный штраф: $\Omega(c) = \|c\|_2^2$ в $\mathbb{R}^{n+1}$.

- Если $\mathcal{F} = C^2[a,b]$ и нежелательна локальная кривизна, то $\hat{f}''$ измеряет изгиб непосредственно. Естественный штраф: $\Omega(\hat{f}) = \|\hat{f}''\|^2_{L^2[a,b]} = \int_a^b (\hat{f}''(x))^2\,dx$.

---

# Регуляризация Тихонова для полиномиальной регрессии

Ищем полином $p(x) = c_0 + c_1 x + \dots + c_n x^n$ степени $n < m$, то есть вектор $c = (c_0,\dots,c_n)^\top \in \mathbb{R}^{n+1}$. Матрица наблюдений $A \in \mathbb{R}^{m\times(n+1)}$, $A_{ij} = x_i^{j-1}$ — матрица Вандермонда, $y \in \mathbb{R}^m$.

Вместо МНК $\min_{c} \|Ac-y\|_2^2$ (число обусловленности $\mathrm{cond}_2(A) = \|A\|_2\|A^*\|_2$ растёт экспоненциально с $n$) решаем задачу Тихонова в $(\mathbb{R}^{n+1}, \|\cdot\|_2)$:

$$\min_{c\,\in\,\mathbb{R}^{n+1}} \Bigl\{ \|Ac - y\|_2^2 + \lambda \|c\|_2^2 \Bigr\}, \quad \lambda > 0$$

**Условия минимума функционала:** $(A^\top A + \lambda I)\,c = A^\top y$; матрица $(A^\top A + \lambda I) \in \mathbb{R}^{(n+1)\times(n+1)}$ невырождена при любом $\lambda > 0$.

---

# Сглаживающий кубический сплайн

**Вариационное свойство.**

Среди всех $\varphi \in C^2[a,b]$ с $\varphi(x_k) = f_k$, $k=0,\dots,n$, естественный кубический сплайн $S$ минимизирует функционал энергии изгиба в $\mathbb{R}$:

$$E(\varphi) = \int_a^b \bigl(\varphi''(x)\bigr)^2\,dx \;\ge\; E(S) \quad \forall\,\varphi$$

---

# Сглаживающий кубический сплайн

Для зашумлённых данных $(x_i, y_i)$, $i=1,\dots,m$, ищем $\hat{f} \in C^2[a,b]$, минимизирующую:

$$\boxed{\mathcal{J}(\hat{f}) = \frac{1}{m}\sum_{i=1}^m \bigl(\hat{f}(x_i) - y_i\bigr)^2 + \lambda \int_a^b \bigl(\hat{f}''(x)\bigr)^2\,dx}$$

Первое слагаемое — среднеквадратичная невязка в $\mathbb{R}$, второе — штраф за кривизну в $L^2[a,b]$.

**Теорема (Reinsch, 1967):** решение — **естественный кубический сплайн** с узлами в точках данных $x_1, \dots, x_m$, не проходящий через данные точно.

---

# Вычисление сглаживающего сплайна

Обозначим $\hat{f} = (\hat{f}(x_1),\dots,\hat{f}(x_m))^\top \in \mathbb{R}^m$ — вектор значений сплайна в узлах, $y = (y_1,\dots,y_m)^\top \in \mathbb{R}^m$.

Из лекций: вторые производные $u_k = S''(x_k)$ связаны со значениями $f_k$ через трёхдиагональную систему $Tu = \rho(f)$, $T \in \mathbb{R}^{(m-1)\times(m-1)}$.

**Матричная форма функционала** (Green & Silverman, 1994): штраф записывается как квадратичная форма $\int_a^b(\hat{f}'')^2\,dx = \hat{f}^\top K \hat{f}$, где $K \in \mathbb{R}^{m\times m}$ — симметричная положительно полуопределённая ленточная матрица, определяемая узлами $\{x_i\}$:

$$\mathcal{J}(\hat{f}) = \frac{1}{m}\|\hat{f} - y\|_2^2 + \lambda\, \hat{f}^\top K \hat{f}$$

**Решение** из условия стационарности $\nabla_{\hat{f}}\,\mathcal{J} = 0$:

$$(I + m\lambda\, K)\,\hat{f} = y, \qquad I,\, K \in \mathbb{R}^{m\times m}$$

- Трехдиагональность $K$ → (метод Томаса) прогонка за $O(m)$ операций
- **Квазилокальность** $(I+m\lambda K)^{-1}$ (какая-то теорема): влияние каждой точки данных убывает экспоненциально с расстоянием

---

# Выбор параметра регуляризации

Оба метода — Тихонов и сглаживающий сплайн — линейны по $y$: $\hat{f}_\lambda = H(\lambda)\,y \in \mathbb{R}^m$, где $H(\lambda) \in \mathbb{R}^{m\times m}$:

$$H(\lambda) = \begin{cases} A(A^\top A + \lambda I)^{-1}A^\top & \text{Тихонов} \\ (I + m\lambda K)^{-1} & \text{сглаживающий сплайн} \end{cases}$$

$H_{ii}(\lambda)$ — диагональные элементы $H(\lambda)$.

## 1. Перекрёстная проверка (LOO-CV)

$$\mathrm{CV}(\lambda) = \frac{1}{m}\sum_{i=1}^m \left(\frac{y_i - \hat{f}_\lambda(x_i)}{1 - H_{ii}(\lambda)}\right)^2 \in \mathbb{R}$$

---

# Выбор параметра регуляризации

## 2. Обобщённая перекрёстная проверка (GCV, Craven & Wahba, 1979)

$$\mathrm{GCV}(\lambda) = \frac{\tfrac{1}{m}\|y - \hat{f}_\lambda\|_2^2}{\left(1 - \tfrac{1}{m}\operatorname{tr} H(\lambda)\right)^2} \in \mathbb{R}$$

Для Тихонова: $\operatorname{tr} H(\lambda) = \sum_{i=1}^r w_i(\lambda) = \sum_{i=1}^r \dfrac{\sigma_i^2}{\sigma_i^2+\lambda}$ — вычисляется через SVD без обращения матриц.

## 3. L-кривая (Hansen, 1992)

| Метод | Оси лог-лог графика |
|---|---|
| Тихонов | $\|Ac_\lambda - y\|_2$ vs. $\|c_\lambda\|_2$ |
| Сглаживающий сплайн | $\|y - \hat{f}_\lambda\|_2$ vs. $\|\hat{f}''_\lambda\|_{L^2[a,b]}$ |

Оптимальное $\lambda^*$ — в точке максимальной кривизны L-образной кривой.

---

# Регуляризация в спектральном базисе

## Регуляризация Тихонова для полиномиальной регрессии

Используя SVD $A = U\Sigma V^\top$ и фильтрующие коэффициенты $w_i(\lambda) = \dfrac{\sigma_i^2}{\sigma_i^2+\lambda} \in [0,1]$ (введены ранее):

$$\hat{c}_\lambda = \sum_{i=1}^r w_i(\lambda)\,\frac{u_i^\top y}{\sigma_i}\,v_i \;\in\; \mathbb{R}^{n+1}$$

## Сглаживающий кубический сплайн

Спектральное разложение $K = V_K D_K V_K^\top \in \mathbb{R}^{m\times m}$, где $D_K = \operatorname{diag}(\mu_1,\dots,\mu_m)$, $\mu_i \ge 0$ — собственные значения $K$, $v_{K,i} \in \mathbb{R}^m$ — ортонормированные собственные векторы. Аналогично вводим $\tilde{w}_i(\lambda) = \dfrac{1}{1+m\lambda\mu_i} \in [0,1]$:

$$(I + m\lambda K)\hat{f} = y \;\;\Rightarrow\;\; \hat{f} = \sum_{i=1}^m \tilde{w}_i(\lambda)\,(v_{K,i}^\top y)\,v_{K,i} \;\in\; \mathbb{R}^m$$

---

# Сравнение подходов

| Критерий | Полином. регрессия + Тихонов | Сглаживающий сплайн |
|---|---|---|
| Класс функций | $\Pi_n \subset C[a,b]$, фиксирован | $C^2[a,b]$, определяется данными |
| Пространство решения | $c \in \mathbb{R}^{n+1}$ | $\hat{f} \in \mathbb{R}^m$ (значения в узлах) |
| Штраф $\Omega$ | $\|c\|_2^2$ в $\mathbb{R}^{n+1}$ | $\int_a^b(\hat{f}'')^2dx$ в $L^2[a,b]$ |
| Параметры | $n$ и $\lambda$ | только $\lambda$ |
| Локальность | глобальная | квазилокальная |
| Сложность | $O(mn^2)$ через SVD | $O(m)$ через метод Томаса |

---

# Заключение

1. **Точная интерполяция в $\bigl(C[a,b],\|\cdot\|_C\bigr)$ неустойчива:**
   $\Lambda_n = \|P_n\|_{C[a,b]\to C[a,b]} \to \infty$ для любой сетки (теорема Бернштейна–Фабера)

2. **Регуляризация:** минимизация $\mathcal{J}(\hat{f}) = \mathcal{D}(\hat{f},y) + \lambda\,\Omega(\hat{f})$; параметр $\lambda>0$ управляет балансом точность - гладкость

3. **Полиномиальная регрессия + Тихонов** в $(\mathbb{R}^{n+1}, \|\cdot\|_2)$: решение $c_\lambda \in \mathbb{R}^{n+1}$ через SVD матрицы Вандермонда $A \in \mathbb{R}^{m\times(n+1)}$; фильтрация через $w_i(\lambda) = \sigma_i^2/(\sigma_i^2+\lambda)$

4. **Сглаживающие сплайны** в $C^2[a,b]$: обобщение вариационного свойства е.с.; решение $\hat{f} \in \mathbb{R}^m$ за $O(m)$ прогонкой; квазилокальность из ленточности $B = I + m\lambda K$

5. **Выбор $\lambda^*$:** LOO-CV, GCV, L-кривая применимы к обоим методам
