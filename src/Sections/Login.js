import React, { Component } from 'react';

import { Navigate } from 'react-router-dom';
import { useState, useEffect } from 'react';

import { variables } from './Variables.js';


export class Login extends Component {

  constructor(props) {
    super(props);

    this.state = {
      token: variables.token,
      email: "",
      password: "",
      error: "",
      isAuthenticated: false,
    }
  }

  componentDidMount() {
    return
  }

  handlePasswordChange = (event) => {
    this.setState({ password: event.target.value });
  }

  handleEmailChange = (event) => {
    this.setState({ email: event.target.value });
  }

  isResponseOk(response) {
    if (response.status >= 200 && response.status <= 299) {
      return response.json();
    } else {
      throw Error(response.statusText);
    }
  }

  login = (event) => {
    event.preventDefault();
    fetch(variables.API_URL + "/accounts/login", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      credentials: "same-origin",
      body: JSON.stringify({ email: this.state.email, password: this.state.password }),
    })
      .then(this.isResponseOk)
      .then((data) => {
        console.log(data);
        this.setState({ token: data.token.key });
        variables.token = data.token.key
        variables.email = data.email
        variables.allow = data.allow
        variables.admin = data.is_admin
      })
      .catch((err) => {
        console.log(err);
        this.setState({ error: "Wrong email or password." });
      });
  }

  logout = () => {
    fetch(variables.API_URL + "/accounts/logout",
      {
        headers: {
          'Content-Type': 'application/json;charset=utf-8',
          'Authorization': `Token ${this.state.token}`,
        },
      }
    )
      .then(this.isResponseOk)
      .then((data) => {
        console.log(data);
        this.setState({ token: null });
        variables.token = null;
      })
      .catch((err) => {
        console.log(err);
      });
  };

  render() {
    if (!this.state.token) {
      return (
        <section className="grid grid-cols-1 lg:grid-cols-2 xl:grid-cols-2 2xl:grid-cols-2 gap-4">
          <div class="absolute inset-0 -z-10 mx-0 max-w-none overflow-hidden">
            <div class="absolute left-1/2 top-0 ml-[-38rem] h-[25rem] w-[81.25rem]">
              <div class="absolute inset-0 bg-gradient-to-r from-[#2590EB] to-[#2563EB] opacity-40 [mask-image:radial-gradient(farthest-side_at_top,white,transparent)]">
                <svg
                  aria-hidden="true"
                  class="absolute inset-x-0 inset-y-[-50%] h-[200%] w-full skew-y-[-18deg] fill-black/40 stroke-black/50 mix-blend-overlay"
                >
                  <defs>
                    <pattern
                      id=":rai:"
                      width="72"
                      height="56"
                      patternUnits="userSpaceOnUse"
                      x="-12"
                      y="4"
                    >
                      <path d="M.5 56V.5H72" fill="none"></path>
                    </pattern>
                  </defs>
                  <rect
                    width="100%"
                    height="100%"
                    stroke-width="0"
                    fill="url(#:rai:)"
                  ></rect>
                  <svg x="-12" y="4" class="overflow-visible">
                    <rect
                      stroke-width="0"
                      width="73"
                      height="57"
                      x="288"
                      y="168"
                    ></rect>
                    <rect
                      stroke-width="0"
                      width="73"
                      height="57"
                      x="144"
                      y="56"
                    ></rect>
                    <rect
                      stroke-width="0"
                      width="73"
                      height="57"
                      x="504"
                      y="168"
                    ></rect>
                    <rect
                      stroke-width="0"
                      width="73"
                      height="57"
                      x="720"
                      y="336"
                    ></rect>
                  </svg>
                </svg>
              </div>
              <svg
                viewBox="0 0 1113 440"
                aria-hidden="true"
                class="absolute left-1/2 top-0 ml-[-19rem] w-[69.5625rem] fill-gray-200 blur-[26px] dark:hidden"
              >
                <path d="M.016 439.5s-9.5-300 434-300S882.516 20 882.516 20V0h230.004v439.5H.016Z"></path>
              </svg>
            </div>
          </div>
          <div className="">
            <div className="flex flex-col items-center lg:w-3/5 justify-center px-6 py-8 mx-auto md:h-screen lg:py-0">
              <a
                href="#"
                className="flex text-center uppercase items-center mb-6 text-2xl font-semibold text-gray-900"
              >
                {/* Логотип
                  <img
                    className="w-8 h-8 mr-2"
                    src="https://flowbite.s3.amazonaws.com/blocks/marketing-ui/logo.svg"
                    alt="logo"
                  />*/}
                Machine learning <br />
                without coding
              </a>
              <ul className="mb-8 space-y-4 text-left text-gray-500">
                <li class="flex items-center space-x-3">
                  <svg
                    class="flex-shrink-0 w-5 h-5 text-green-500ы0"
                    fill="currentColor"
                    viewBox="0 0 20 20"
                    xmlns="http://www.w3.org/2000/svg"
                  >
                    <path
                      fill-rule="evenodd"
                      d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z"
                      clip-rule="evenodd"
                    ></path>
                  </svg>
                  <span>
                    Распределённая архитектура: система может быть запущена на
                    вычислительной системе любой сложности.
                  </span>
                </li>
                <li class="flex items-center space-x-3">
                  <svg
                    class="flex-shrink-0 w-5 h-5 text-green-500"
                    fill="currentColor"
                    viewBox="0 0 20 20"
                    xmlns="http://www.w3.org/2000/svg"
                  >
                    <path
                      fill-rule="evenodd"
                      d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z"
                      clip-rule="evenodd"
                    ></path>
                  </svg>
                  <span>
                    Конфигурируемость: система позволяет настраивать параметры
                    алгоритмов для подходящих конкретным задачам.
                  </span>
                </li>
                <li class="flex items-center space-x-3">
                  <svg
                    class="flex-shrink-0 w-5 h-5 text-green-500"
                    fill="currentColor"
                    viewBox="0 0 20 20"
                    xmlns="http://www.w3.org/2000/svg"
                  >
                    <path
                      fill-rule="evenodd"
                      d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z"
                      clip-rule="evenodd"
                    ></path>
                  </svg>
                  <span>
                    Возможность работать с крупными объёмами данных, а также
                    высокая производительность и точность получаемых результатов
                  </span>
                </li>
              </ul>
            </div>
          </div>
          <div className="">
            <div className="flex items-center justify-center px-3 py-8 mx-auto md:h-screen lg:py-0">
              <div className="w-full bg-white rounded-lg drop-shadow-md dark:border md:mt-0 sm:max-w-md xl:p-0">
                <div className="p-6 space-y-4 md:space-y-6 sm:p-8">
                  <h1 className="text-xl font-bold leading-tight tracking-tight text-gray-900 md:text-2xl">
                    Войти в личный кабинет
                  </h1>
                  {/* Авторизация пользователя */}
                  <form className="space-y-4 md:space-y-6" onSubmit={this.login}>
                    <div>
                      <label
                        for="email"
                        className="block mb-2 text-sm font-medium text-gray-900"
                      >
                        Email
                      </label>
                      <input
                        type="email"
                        id="email"
                        name="email"
                        className="bg-gray-50 border border-gray-300 text-gray-900 sm:text-sm rounded-lg focus:ring-primary-600 focus:border-primary-600 block w-full p-2.5"
                        value={this.state.email} onChange={this.handleEmailChange}
                        placeholder="Адрес электронной почты"
                      />
                    </div>
                    <div>
                      <label
                        for="password"
                        className="block mb-2 text-sm font-medium text-gray-900"
                      >
                        Пароль
                      </label>
                      <input
                        type="password"
                        id="password"
                        name="password"
                        value={this.state.password} onChange={this.handlePasswordChange}
                        placeholder="••••••••"
                        className="bg-gray-50 border border-gray-300 text-gray-900 sm:text-sm rounded-lg focus:ring-primary-600 focus:border-primary-600 block w-full p-2.5 "
                      />
                    </div>
                    {this.state.error &&
                      <small className="text-danger">
                        {this.state.error}
                      </small>
                    }
                    <button type="submit" className="w-full text-white bg-blue-600 hover:bg-blue-700 focus:ring-4 focus:outline-none focus:ring-blue-300 font-medium rounded-lg text-sm px-5 py-2.5 text-center"
                    >Войти</button>
                  </form>
                </div>
              </div>
            </div>
          </div>
        </section>
      );
    } else if (variables.allow === 1){
      return <Navigate push to="/ddi_review" />
    } else {
      return <Navigate push to="/tematic_review" />
    }
  }
}

