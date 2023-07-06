import './App.css';
import { useEffect } from 'react';

// Импорт страниц
import { Main } from './Sections/Main.js';
import { Login } from './Sections/Login.js';
import { TematicReview } from './Sections/TematicReview.js';

import { BrowserRouter as Router, Route, Routes } from 'react-router-dom';

function App() {
  useEffect(() => {
    // 👇️ adding multiple classes to the body element
    document.body.classList.add(
      'bg-gray-100',
      'dark:bg-gray-900',
    );
  }, []);

  return (
    <Router>
      <Routes>
        <Route path='/' element={<Main />} />
        <Route path='/login' element={<Login />} />
        <Route path='/tematic_review' element={<TematicReview />} />
      </Routes>
    </Router>
  );

}

export default App;
