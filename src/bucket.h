// A simple integer‐score bucket for peeling.
// extractMin=true ⇒ pop lowest‐score first (for tip‐peeling)
// extractMin=false⇒ pop highest‐score first (for complement‐degeneracy)
struct Bucket {
    int n, maxScore;
    bool extractMin;
    int cur;               // current bucket index we’re scanning
    vector<int> head;      // head[s] = first item in bucket s, or -1
    vector<int> next;      // next[i] = next item after i in its bucket, or -1
    vector<int> score;     // score[i] = current score of item i
  
    // ctor: n items, scores in [0..maxS], extractMin mode
    Bucket(int n_, int maxS, bool extractMin_)
      : n(n_), maxScore(maxS), extractMin(extractMin_),
        head(maxS+1, -1), next(n_,-1), score(n_,0),
        cur(extractMin_?0:maxS)
    {}
  
    // Build the initial buckets from initScore[0..n)
    void initialize(const vector<int>& initScore) {
      for (int i = 0; i < n; i++) score[i] = initScore[i];
      for (int s = 0; s <= maxScore; s++) head[s] = -1;
      for (int i = 0; i < n; i++) {
        int s = score[i];
        next[i]   = head[s];
        head[s]   = i;
      }
      cur = extractMin ? 0 : maxScore;
    }
  
    // Pop next valid item, or -1 if empty
    int pop() {
      if (extractMin) {
        while (cur<=maxScore && head[cur]==-1) cur++;
        if (cur>maxScore) return -1;
      } else {
        while (cur>=0       && head[cur]==-1) cur--;
        if (cur<0) return -1;
      }
      int u = head[cur];
      head[cur] = next[u];
      if (score[u] != cur) return pop();   // stale, skip
      return u;
    }
  
    // Change item i’s score by +delta or –delta (delta = ±1)
    void update(int i, int delta) {
      int old = score[i];
      int now = old + delta;
      score[i] = now;
      next[i]  = head[now];
      head[now] = i;
      if (extractMin) cur = min(cur, now);
      else           cur = max(cur, now);
    }
  
    int getScore(int i) const { return score[i]; }
  };
  